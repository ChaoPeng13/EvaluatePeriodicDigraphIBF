/**
*
*/
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <future>
#include <chrono>
#include <Psapi.h>

#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>
#pragma comment(lib,"E:\\googletest-master\\googletest\\msvc\\gtest\\Debug\\gtestd.lib")
#endif

#ifdef VLDTEST
#include <vld.h>
#endif

#include "SchedulabilityAnalysis.h"
#include "RandomGenerator.h"
#include "FileWriter.h"
#include "FileReader.h"
#include "StateflowExample.h"
#include "ResponseTimeAnalysis.h"
#include "Timer.h"
#include "GraphAlgorithms.h"
#include "Utility.h"
#include "MaxPlusAlgebra.h"
#include "DigraphExperiments.h"
#include "FSMExperiments.h"
#include "DigraphExample.h"

//#pragma comment(linker, "/STACK:4000000")
//#pragma comment(linker, "/HEAP:2000000000")

// To ensure correct resolution of symbols, add Psapi.lib to TARGETLIBS
// and compile with -DPSAPI_VERSION=1

void PrintMemoryInfo( DWORD processID )
{
	HANDLE hProcess;
	PROCESS_MEMORY_COUNTERS pmc;

	// Print the process identifier.

	printf( "\nProcess ID: %u\n", processID );

	// Print information about the memory usage of the process.

	hProcess = OpenProcess(  PROCESS_QUERY_INFORMATION |
		PROCESS_VM_READ,
		FALSE, processID );
	if (NULL == hProcess)
		return;

	if ( GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc)) )
	{
		printf( "\tPageFaultCount: 0x%08X\n", pmc.PageFaultCount );
		printf( "\tPeakWorkingSetSize: 0x%08X\n", 
			pmc.PeakWorkingSetSize );
		printf( "\tWorkingSetSize: 0x%08X\n", pmc.WorkingSetSize );
		printf( "\tQuotaPeakPagedPoolUsage: 0x%08X\n", 
			pmc.QuotaPeakPagedPoolUsage );
		printf( "\tQuotaPagedPoolUsage: 0x%08X\n", 
			pmc.QuotaPagedPoolUsage );
		printf( "\tQuotaPeakNonPagedPoolUsage: 0x%08X\n", 
			pmc.QuotaPeakNonPagedPoolUsage );
		printf( "\tQuotaNonPagedPoolUsage: 0x%08X\n", 
			pmc.QuotaNonPagedPoolUsage );
		printf( "\tPagefileUsage: 0x%08X\n", pmc.PagefileUsage ); 
		printf( "\tPeakPagefileUsage: 0x%08X\n", 
			pmc.PeakPagefileUsage );
	}

	CloseHandle( hProcess );
}

void testGetHeapSize() {
	// Get the list of process identifiers.

	DWORD aProcesses[1024], cbNeeded, cProcesses;
	unsigned int i;

	if ( !EnumProcesses( aProcesses, sizeof(aProcesses), &cbNeeded ) )
	{
		cerr << "Error Here!" << endl;
		exit(EXIT_FAILURE);
	}

	// Calculate how many process identifiers were returned.

	cProcesses = cbNeeded / sizeof(DWORD);

	// Print the memory usage for each process

	for ( i = 0; i < cProcesses; i++ )
	{
		PrintMemoryInfo( aProcesses[i] );
	}

}

bool usingTooMuchMemory()
{
	PROCESS_MEMORY_COUNTERS pmc;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) {
		cout << pmc.WorkingSetSize << endl;
		return pmc.WorkingSetSize > 0x80000000u; // 2GB working set
	}
	return false;
}

void testMaxHeapMemory() {
	int iter = 0;
	try {
		while(true) {
			int* ptr = new int[1024*1024];
			cout << ++iter << endl;
		}
		cout << usingTooMuchMemory() << endl;
	}
	catch(bad_alloc& ba)
	{
		cout << "sizeof(int) = " << sizeof(int) << endl;
		cout<<iter*4/1024<<endl;
	}
}

void testRBFIBF(string directory, string file, int minRun, int maxRun, int stepRun) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int exception = 0;
	Timer timer;
	Digraph* digraph;
	double tibf0 = 0;
	double tibf1 = 0;
	double tibf2 = 0;
	double tibf3 = 0;
	double trbf0 = 0;
	double trbf1 = 0;
	double trbf2 = 0;
	double trbf3 = 0;

	try {

		for (int run=minRun; run<=maxRun; run+=stepRun) {
			/*
			string name = directory + "\\digraph" + Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			FileReader::DotFileReader(digraph,1,p);
			cout << name << endl;
			*/

			digraph = DigraphExample::generateDigraph6();

			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			/*
			if (!digraph->strongly_connected) {
			cerr << "The generated digraph has to be strongly connected." << endl;
			exit(EXIT_FAILURE);
			}
			*/

			digraph->calculate_period_gcd();
			digraph->calculate_all_gcd();
			digraph->calculate_linear_factor();
			//digraph->tf = 50000; // tf is large enough to calculate linear defect

			// calculate the linear period
			digraph->unit_digraph = new UnitDigraph(digraph);
			digraph->unit_digraph->prepare3(false);


			cout << "lfac=" << digraph->linear_factor << "\tlper=" << digraph->unit_digraph->lper << endl;

			int tf = 100000;

			timer.start();
			digraph->calculate_ibf_without_periodicity(tf,digraph->ibf_map2);
			timer.end();
			tibf0 += timer.getTime();

			timer.start();
			digraph->calculate_ibf_without_periodicity_fast(tf,digraph->ibf_map2_fast);
			timer.end();
			tibf1 += timer.getTime();

			timer.start();
			digraph->calculate_ibf_without_periodicity_DP(tf,digraph->ibf_map2_DP);
			timer.end();
			tibf2 += timer.getTime();

			map<int,int> ibf_map3;
			timer.start();
			digraph->calculate_ibf_with_periodicity_DP(tf,ibf_map3);
			timer.end();
			tibf3 += timer.getTime();

			timer.start();
			digraph->calculate_rbf_without_periodicity(tf,digraph->rbf_map2);
			timer.end();
			trbf0 += timer.getTime();

			timer.start();
			digraph->calculate_rbf_without_periodicity_fast(tf,digraph->rbf_map2_fast);
			timer.end();
			trbf1 += timer.getTime();

			timer.start();
			digraph->calculate_rbf_without_periodicity_DP(tf,digraph->rbf_map2_DP);
			timer.end();
			trbf2 += timer.getTime();

			map<int,int> rbf_map3;
			timer.start();
			digraph->calculate_rbf_with_periodicity_DP(tf,rbf_map3);
			timer.end();
			trbf3 += timer.getTime();

			map<int,int> dbf_map2;
			digraph->calculate_dbf_without_periodicity_DP(tf,dbf_map2);

			cout << "Calculating time for ibf..." << endl;
			cout << "tibf0=" << tibf0 << "\ttibf1=" << tibf1 << "\ttibf2=" << tibf2 << "\ttibf3=" << tibf3 << endl;
			cout << "trbf0=" << trbf0 << "\ttrbf1=" << trbf1 << "\ttrbf2=" << trbf2 << "\ttrbf3=" << trbf3 << endl;

			if (true) {
				/*
				for (int t=0; t<=100; t+= digraph->aGCD) {
				int ibf0 = digraph->ibf2(t);
				int ibf1 = digraph->ibf2_fast(t);
				cout << "ibf(" << t << ")="<< ibf0 << "\t" << ibf1 << endl;


				if (ibf0 != ibf1) {
				cerr << "At time " << t << " the two algorithms have different values: " << ibf0 << "\t" << ibf1 << endl;
				exit(EXIT_FAILURE);
				}

				}
				*/

				for(map<int,int>::iterator iter = digraph->ibf_map2_fast.begin(); iter != digraph->ibf_map2_fast.end(); iter++) {
					if (iter->first <= tf) {
						int t = iter->first;
						int ibf0 = digraph->ibf2(t);
						int ibf1 = digraph->ibf2_fast(t);
						int ibf2 = digraph->ibf2(t,digraph->ibf_map2_DP);
						int ibf3 = digraph->ibf2(t,ibf_map3);
						//cout << "ibf(" << t << ")="<< ibf0 << "\t" << ibf2 << endl;

						if (ibf0 != ibf1) {
							cerr << "(1.ibf) At time " << t << " the two algorithms have different values: " << ibf0 << "\t" << ibf1 << endl;
							exit(EXIT_FAILURE);
						}

						if (ibf0 != ibf2) {
							cerr << "(2.ibf) At time " << t << " the two algorithms have different values: " << ibf0 << "\t" << ibf2 << endl;
							exit(EXIT_FAILURE);
						}

						if (ibf0 != ibf3) {
							cerr << "(3.ibf) At time " << t << " the two algorithms have different values: " << ibf0 << "\t" << ibf3 << endl;
							exit(EXIT_FAILURE);
						}
					}
				}

				/*
				for(map<int,int>::iterator iter = digraph->rbf_map2.begin(); iter != digraph->rbf_map2.end(); iter++) {
				if (iter->first <= 58) {
				cout << "rbf(" << iter->first << ")="<< iter->second << endl;
				}
				}
				*/

				for(map<int,int>::iterator iter = digraph->rbf_map2_fast.begin(); iter != digraph->rbf_map2_fast.end(); iter++) {
					if (iter->first <= tf-100) {
						int t = iter->first;
						int rbf0 = digraph->rbf2(t);
						int rbf1 = digraph->rbf2_fast(t);
						int rbf2 = digraph->rbf2(t,digraph->rbf_map2_DP);
						int rbf3 = digraph->rbf2(t,rbf_map3);
						//cout << "rbf(" << t << ")="<< rbf0 << "\t" << rbf2 << endl;

						if (rbf0 != rbf1) {
							cerr << "(1) At time " << t << " the two algorithms have different values: " << rbf0 << "\t" << rbf1 << endl;
							exit(EXIT_FAILURE);
						}

						if (rbf0 != rbf2) {
							cerr << "(2) At time " << t << " the two algorithms have different values: " << rbf0 << "\t" << rbf2 << endl;
							exit(EXIT_FAILURE);
						}

						if (rbf0 != rbf3) {
							cerr << "(3) At time " << t << " the two algorithms have different values: " << rbf0 << "\t" << rbf3 << endl;
							exit(EXIT_FAILURE);
						}
					}
				}

				for(map<int,int>::iterator iter = dbf_map2.begin(); iter != dbf_map2.end(); iter++) {
					if (iter->first <= tf-100) {
						int t = iter->first;
						int dbf0 = digraph->calculate_dbf(t);
						int dbf1 = digraph->dbf2(t,dbf_map2);
						cout << "dbf(" << t << ")="<< dbf0 << "\t" << dbf1 << endl;

						if (dbf0 != dbf1) {
							cerr << "(1) At time " << t << " the two algorithms have different values: " << dbf0 << "\t" << dbf1 << endl;
							exit(EXIT_FAILURE);
						}
					}
				}
			}
		}
	}
	catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
		delete digraph;
	}

	cout << "Calculating time for ibf..." << endl;
	cout << "tibf0=" << tibf0 << "\ttibf1=" << tibf1 << "\ttibf2=" << tibf2 << "\ttibf3=" << tibf3 << endl;
	cout << "trbf0=" << trbf0 << "\ttrbf1=" << trbf1 << "\ttrbf2=" << trbf2 << "\ttrbf3=" << trbf3 << endl;

	fout.close();
}

int main(int argc, char* argv[])
{
#ifdef GOOGLETEST
	// Run Google test
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
#endif
	string index = "20161220";
	string result = "Results\\";
	string suffix = ".res";
	//int period_choice = RandomGenerator::ExactAnalysis3;
	//int period_choice = RandomGenerator::ExactAnalysis4;
	//int period_choice = RandomGenerator::ExactAnalysis5;
	//int period_choice = RandomGenerator::ExactAnalysis6;
	//int period_choice = RandomGenerator::ExactAnalysis7;
	//int period_choice = RandomGenerator::ExactAnalysis8;
	//int period_choice = RandomGenerator::ExactAnalysis9;
	//int period_choice = RandomGenerator::ApproximateAnalysis5;
	//int period_choice = RandomGenerator::ApproximateAnalysis9;
	//int period_choice = RandomGenerator::ApproximateAnalysis10;
	//int period_choice = RandomGenerator::ApproximateAnalysis11;
	//int period_choice = RandomGenerator::ApproximateAnalysis12;
	//int period_choice = RandomGenerator::ApproximateAnalysis13;

	//int period_choice = RandomGenerator::ApproximateAnalysis0;
	//int period_choice = RandomGenerator::ApproximateAnalysis6;
	//int period_choice = RandomGenerator::ApproximateAnalysis7;

	// myProperty = true, do exact analysis
	// otherwise, do approximate analysis
	bool myProperty = false; 
	int times = 10;
	int period_choice = -1;

	//testGetHeapSize();

	// MethodChoice mc = EXACTSO;
	// MethodChoice mc = EXACTAO1;
	// MethodChoice mc = EXACTAO2;
	//MethodChoice mc = IBFSO;
	// MethodChoice mc = RBFSO;
	// MethodChoice mc = IBFSOP;
	// MethodChoice mc = RBFSOP;
	// MethodChoice mc = RBFAOFSM;
	// MethodChoice mc = IBFAOFSM;
	// MethodChoice mc = RBFAOFSMDP;
	// MethodChoice mc = IBFAOFSMDP;
	// MethodChoice mc = IBFAO1G;
	// MethodChoice mc = RBFAO1G;
	// MethodChoice mc = IBFAO1DP;
	// MethodChoice mc = RBFAO1DP;
	// MethodChoice mc = RBFAO1L;
	// MethodChoice mc = IBFAO2G;
	// MethodChoice mc = RBFAO2G;
	// MethodChoice mc = IBFAO2DP;
	//MethodChoice mc = RBFAO2DP;
	// MethodChoice mc = RBFAO2L;
	// MethodChoice mc = LUBIBF;
	// MethodChoice mc = LUBRBF;
	// MethodChoice mc = LUBCSUM;

	//for (int mc = IBFSO; mc <= LUBCSUM; mc++) {

	int uIndex = 5;
	int nIndex = 2;
	int sIndex = 5;

	//for (int iii = 0; iii <8 ; iii ++)
	//cout << RandomGenerator::STATEFLOW_PERIOD5[iii] << endl;

	if (false) {
		period_choice = RandomGenerator::ExactAnalysis8;
		srand (time(NULL));
		
		if (argc != 9) {
			cerr << "Error parameter number => " << argc << endl;
			cerr << "Info: *.exe 20 20 5 60 0 999 Util All" << endl;
			cerr << "      *.exe 5 26 30 30 0 999 Num All" << endl;
			exit(EXIT_FAILURE);
		}

		int minTaskNum = atoi(argv[1]);
		int maxTaskNum = atoi(argv[2]);
		int minUtil = atoi(argv[3]);
		int maxUtil = atoi(argv[4]);
		int minRun = atoi(argv[5]);
		int maxRun = atoi(argv[6]); 
		string testName(argv[7]);
		string fileName(argv[8]); 

		if (testName == "Util") {
			// varying util
			FSMExperiments::generateRandomSystemsForExactAnalysis("ExactAnalysisUtil"+index+"Period"+Utility::int_to_string(period_choice),period_choice,minTaskNum,maxTaskNum,1,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			FSMExperiments::SchedAnalysis(EXACTSO,result+"ExactAnalysisUtil"+index+"Period"+Utility::int_to_string(period_choice)+"file"+fileName,"ExactAnalysisUtil"+index+"Period"+Utility::int_to_string(period_choice),minTaskNum,maxTaskNum,1,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSO,result+"ExactAnalysisUtil"+index+"Period"+Utility::int_to_string(period_choice)+"file"+fileName,"ExactAnalysisUtil"+index+"Period"+Utility::int_to_string(period_choice),minTaskNum,maxTaskNum,1,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSOP,result+"ExactAnalysisUtil"+index+"Period"+Utility::int_to_string(period_choice)+"file"+fileName,"ExactAnalysisUtil"+index+"Period"+Utility::int_to_string(period_choice),minTaskNum,maxTaskNum,1,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
		} else {
			// varying task number
			FSMExperiments::generateRandomSystemsForExactAnalysis("ExactAnalysisNum"+index+"Period"+Utility::int_to_string(period_choice),period_choice,minTaskNum,maxTaskNum,3,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			FSMExperiments::SchedAnalysis(EXACTSO,result+"ExactAnalysisNum"+index+"Period"+Utility::int_to_string(period_choice)+"file"+fileName,"ExactAnalysisNum"+index+"Period"+Utility::int_to_string(period_choice),minTaskNum,maxTaskNum,3,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSO,result+"ExactAnalysisNum"+index+"Period"+Utility::int_to_string(period_choice)+"file"+fileName,"ExactAnalysisNum"+index+"Period"+Utility::int_to_string(period_choice),minTaskNum,maxTaskNum,3,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSOP,result+"ExactAnalysisNum"+index+"Period"+Utility::int_to_string(period_choice)+"file"+fileName,"ExactAnalysisNum"+index+"Period"+Utility::int_to_string(period_choice),minTaskNum,maxTaskNum,3,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
		}
	}

	if (false) {
		//period_choice = RandomGenerator::ExactAnalysis8;
		srand (time(NULL));
		
		if (argc != 9) {
			cerr << "Error parameter number => " << argc << endl;
			cerr << "Info: *.exe 20 20 5 60 0 999 Util All" << endl;
			cerr << "      *.exe 5 60 30 30 0 999 Num All" << endl;
			exit(EXIT_FAILURE);
		}

		int minTaskNum = atoi(argv[1]);
		int maxTaskNum = atoi(argv[2]);
		int minUtil = atoi(argv[3]);
		int maxUtil = atoi(argv[4]);
		int minRun = atoi(argv[5]);
		int maxRun = atoi(argv[6]); 
		string testName(argv[7]);
		string fileName(argv[8]); 

		period_choice = RandomGenerator::ApproximateAnalysis13;
		// evaluate the periodicity property
		if (testName == "Util") {
			FSMExperiments::generateRandomSystemsForApproximateAnalysisWithPeriodicity2("PApproximateAnalysisUtil"+index,period_choice,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			// FSMExperiments::SchedAnalysis(mc,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			FSMExperiments::SchedAnalysis(LUBIBF,result+"PApproximateAnalysisUtil"+index+"file"+fileName,"PApproximateAnalysisUtil"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBRBF,result+"PApproximateAnalysisUtil"+index+"file"+fileName,"PApproximateAnalysisUtil"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBCSUM,result+"PApproximateAnalysisUtil"+index+"file"+fileName,"PApproximateAnalysisUtil"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSO,result+"PApproximateAnalysisUtil"+index+"file"+fileName,"PApproximateAnalysisUtil"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSOP,result+"PApproximateAnalysisUtil"+index+"file"+fileName,"PApproximateAnalysisUtil"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSO,result+"PApproximateAnalysisUtil"+index+"file"+fileName,"PApproximateAnalysisUtil"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSOP,result+"PApproximateAnalysisUtil"+index+"file"+fileName,"PApproximateAnalysisUtil"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
		} else {
			FSMExperiments::generateRandomSystemsForApproximateAnalysisWithPeriodicity2("PApproximateAnalysisNum"+index,period_choice,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			//FSMExperiments::SchedAnalysis(mc,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,30,30,5,90,90,10,0,times-1);
			FSMExperiments::SchedAnalysis(LUBIBF,result+"PApproximateAnalysisNum"+index+"file"+fileName,"PApproximateAnalysisNum"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBRBF,result+"PApproximateAnalysisNum"+index+"file"+fileName,"PApproximateAnalysisNum"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBCSUM,result+"PApproximateAnalysisNum"+index+"file"+fileName,"PApproximateAnalysisNum"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSO,result+"PApproximateAnalysisNum"+index+"file"+fileName,"PApproximateAnalysisNum"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSOP,result+"PApproximateAnalysisNum"+index+"file"+fileName,"PApproximateAnalysisNum"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSO,result+"PApproximateAnalysisNum"+index+"file"+fileName,"PApproximateAnalysisNum"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSOP,result+"PApproximateAnalysisNum"+index+"file"+fileName,"PApproximateAnalysisNum"+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,90,90,10,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
		}
	}

	
	// varying util
	//FSMExperiments::SchedAnalysis(IBFAOFSMDP,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,15,15,5,90,90,10,0,0);

	if (false) {
		period_choice = RandomGenerator::ApproximateAnalysis5;
		if (true) {
			//FSMExperiments::generateRandomSystemsForApproximateAnalysis("ApproximateAnalysisUtil"+index,period_choice,20,20,1,5,60,5,90,90,10,times);
			FSMExperiments::SchedAnalysis(RBFSO,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			//FSMExperiments::SchedAnalysis(mc,result+"ApproximateAnalysisUtil2"+index,"ApproximateAnalysisUtil"+index,20,20,1,30,30,5,90,90,10,153,153);
			//FSMExperiments::generateRandomSystemsForApproximateAnalysis("ApproximateAnalysisUtil"+index,period_choice,2,2,1,5,5,5,90,90,10,1);
			//FSMExperiments::SchedAnalysis(mc,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,2,2,1,5,5,5,90,90,10,0,0);
			//FSMExperiments::SchedAnalysis(mc,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,25,25,5,90,90,10,560,560);
			//int FSMUtil = 5;
			//FSMExperiments::SchedAnalysis(mc,result+"ApproximateAnalysisUtil"+Utility::int_to_string(FSMUtil)+index,"ApproximateAnalysisUtil"+index,20,20,1,FSMUtil,FSMUtil+5,5,90,90,10,0,times-1);
			/*
			FSMExperiments::SchedAnalysis(RBFAOFSMDP,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFAO1DP,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO1DP,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFAO2DP,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO2DP,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSO,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSO,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFAOFSM,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAOFSM,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFAOFSMDP,result+"ApproximateAnalysisUtil"+index,"ApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			*/

		} else {
			// varying task number
			FSMExperiments::generateRandomSystemsForApproximateAnalysis("ApproximateAnalysisNum"+index,period_choice,5,60,5,40,40,5,90,90,10,times);
			//FSMExperiments::SchedAnalysis(mc,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			//int FSMNum = 5;
			//FSMExperiments::SchedAnalysis(mc,result+"ApproximateAnalysisNum"+Utility::int_to_string(FSMNum)+index,"ApproximateAnalysisNum"+index,FSMNum,FSMNum+5,5,40,40,5,90,90,10,0,times-1);                                                                                                              
			FSMExperiments::SchedAnalysis(RBFAOFSMDP,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(IBFAO1DP,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);  
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(RBFAO1DP,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);  
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(IBFAO2DP,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1); 
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(RBFAO2DP,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);  
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(IBFSO,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);     
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(RBFSO,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);     
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(IBFAOFSM,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);  
			SchedulabilityAnalysis::clear_all();                                                                                                               
			FSMExperiments::SchedAnalysis(RBFAOFSM,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);  
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFAOFSMDP,result+"ApproximateAnalysisNum"+index,"ApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);  
			SchedulabilityAnalysis::clear_all(); 
		}
	}

	int choice = 2;
	if (false) {
		srand (time(NULL));
		
		if (argc != 11) {
			cerr << "Error parameter number => " << argc << endl;
			cerr << "Info: *.exe 20 20 40 40 25 500 0 9 Scalability All" << endl;
			exit(EXIT_FAILURE);
		}

		int minTaskNum = atoi(argv[1]);
		int maxTaskNum = atoi(argv[2]);
		int minUtil = atoi(argv[3]);
		int maxUtil = atoi(argv[4]);
		int minState = atoi(argv[5]);
		int maxState = atoi(argv[6]);
		int minRun = atoi(argv[7]);
		int maxRun = atoi(argv[8]); 
		string testName(argv[9]);
		string fileName(argv[10]); 

		period_choice = RandomGenerator::ApproximateAnalysis13;

		FSMExperiments::generateRandomSystemsForScalabilityTest("ApproximateAnalysisState"+index,period_choice,minTaskNum,maxTaskNum,1,minUtil,maxUtil,5,90,90,10,minState,maxState,25, minRun, maxRun);
		FSMExperiments::SchedAnalysisForScalabilityTest(RBFSO,myProperty,result+"ApproximateAnalysisState"+index+"file"+fileName,"ApproximateAnalysisState"+index, choice, minTaskNum, minUtil, 90, minState, maxState, 25, minRun, maxRun);
		SchedulabilityAnalysis::clear_all();
		FSMExperiments::SchedAnalysisForScalabilityTest(IBFSO,myProperty,result+"ApproximateAnalysisState"+index+"file"+fileName,"ApproximateAnalysisState"+index, choice, minTaskNum, minUtil, 90, minState, maxState, 25, minRun, maxRun);
		SchedulabilityAnalysis::clear_all();
		FSMExperiments::SchedAnalysisForScalabilityTest(LUBIBF,myProperty,result+"ApproximateAnalysisState"+index+"file"+fileName,"ApproximateAnalysisState"+index, choice, minTaskNum, minUtil, 90, minState, maxState, 25, minRun, maxRun);
		SchedulabilityAnalysis::clear_all();
		FSMExperiments::SchedAnalysisForScalabilityTest(LUBRBF,myProperty,result+"ApproximateAnalysisState"+index+"file"+fileName,"ApproximateAnalysisState"+index, choice, minTaskNum, minUtil, 90, minState, maxState, 25, minRun, maxRun);
		SchedulabilityAnalysis::clear_all();
		//FSMExperiments::SchedAnalysisForScalabilityTest(LUBCSUM,myProperty,result+"ApproximateAnalysisState"+index,"ApproximateAnalysisState"+index, choice, 20, 40, 90, 25, 500, 25, 0, times-1);
		//SchedulabilityAnalysis::clear_all();
	}

	if (false) {
		period_choice = RandomGenerator::ApproximateAnalysis13;
		index += "SN15#2#";
		FSMExperiments::generateRandomSystemsForWeightedSchedAnalysis("WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
		FSMExperiments::WeightedSchedulabilityAnalysis(IBFSO,"WeightedSched"+index, result+"WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
		SchedulabilityAnalysis::clear_all();
		FSMExperiments::WeightedSchedulabilityAnalysis(RBFSO,"WeightedSched"+index, result+"WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
		SchedulabilityAnalysis::clear_all();
		FSMExperiments::WeightedSchedulabilityAnalysis(LUBIBF,"WeightedSched"+index, result+"WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
		SchedulabilityAnalysis::clear_all();
		FSMExperiments::WeightedSchedulabilityAnalysis(LUBRBF,"WeightedSched"+index, result+"WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
		SchedulabilityAnalysis::clear_all();
		FSMExperiments::WeightedSchedulabilityAnalysis(LUBCSUM,"WeightedSched"+index, result+"WeightedSched"+index,period_choice,5,60,5,0,100,90,90,10,times);
		SchedulabilityAnalysis::clear_all();
	}
	//FSMExperiments::SchedAnalysis(RBFAO1LSCC,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
	//SchedulabilityAnalysis::clear_all();

	if (false) {
		period_choice = RandomGenerator::ApproximateAnalysis13;
		index += "SN15#2#";
		// evaluate the periodicity property
		if (true) {
			FSMExperiments::generateRandomSystemsForApproximateAnalysisWithPeriodicity("PApproximateAnalysisUtil"+index,period_choice,20,20,1,5,60,5,90,90,10,times);
			// FSMExperiments::SchedAnalysis(mc,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			FSMExperiments::SchedAnalysis(LUBIBF,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBRBF,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBCSUM,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSO,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSOP,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO1DP,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO1L,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			//FSMExperiments::SchedAnalysis(RBFAO1LSCC,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			//SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO2DP,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO2L,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSO,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSOP,result+"PApproximateAnalysisUtil"+index,"PApproximateAnalysisUtil"+index,20,20,1,5,60,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
		} else {
			FSMExperiments::generateRandomSystemsForApproximateAnalysisWithPeriodicity("PApproximateAnalysisNum"+index,period_choice,5,60,5,40,40,5,90,90,10,times);
			//FSMExperiments::SchedAnalysis(mc,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,30,30,5,90,90,10,0,times-1);
			FSMExperiments::SchedAnalysis(LUBIBF,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBRBF,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(LUBCSUM,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFSOP,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO1DP,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO1L,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			//FSMExperiments::SchedAnalysis(RBFAO1LSCC,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			//SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO2DP,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(RBFAO2L,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSO,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
			FSMExperiments::SchedAnalysis(IBFSOP,result+"PApproximateAnalysisNum"+index,"PApproximateAnalysisNum"+index,5,60,5,40,40,5,90,90,10,0,times-1);
			SchedulabilityAnalysis::clear_all();
		}
	}

	/*
	if (true) {
	DigraphExperiments::generateDigraphs("Digraph"+index,20,20,5,70,70,5,1);
	DigraphExperiments::testDigraphs("Digraph"+index,result+"DigraphUtil"+index,20,20,5,70,70,5,1);
	}
	else {
	DigraphExperiments::generateDigraphs("Digraph"+index,5,40,5,70,70,5,10);
	DigraphExperiments::testDigraphs("Digraph"+index,result+"DigraphNum"+index,5,40,5,70,70,5,10);

	DigraphExperiments::generateDigraphs("Digraph"+index,5,40,5,80,80,5,10);
	DigraphExperiments::testDigraphs("Digraph"+index,result+"DigraphNum"+index,5,40,5,80,80,5,10);

	DigraphExperiments::generateDigraphs("Digraph"+index,5,40,5,90,90,5,10);
	DigraphExperiments::testDigraphs("Digraph"+index,result+"DigraphNum"+index,5,40,5,90,90,5,10);
	}
	*/

	// test 
	
#if 0
	times = 1000;
	for (int numNode = 5; numNode <= 25; numNode += 5) {
		DigraphExperiments::generateOneDigraphForTestingImprovemntOnIBFCalculation2("Digraphs"+index+"Node"+Utility::int_to_string(numNode),50,numNode,times);
		DigraphExperiments::testImprovemntOnIBFCalculationInOneDigrah2("Digraphs"+index+"Node"+Utility::int_to_string(numNode),result+"DigraphTimeCalIBF"+index+"Node"+Utility::int_to_string(numNode),50,numNode,0,times-1,1);
	}
#endif

#if 1
	times = 1000;
	DigraphExperiments::generateOneDigraphForTestingImprovemntOnIBFCalculation2("Digraphs"+index+"NodeRandom",50,10,25,times);
	DigraphExperiments::testImprovemntOnIBFCalculationInOneDigrah2("Digraphs"+index+"NodeRandom",result+"DigraphTimeCalIBF"+index+"NodeRandom",50,1025,0,times-1,1);
#endif

#if 0
	times = 1000;
	bool deadlineProperty = true;
	//srand (time(NULL));

	//DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2("ImprovementUtil"+index,deadlineProperty, 10, 25,20,20,1,50,95,5,0,times-1);
	
	DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"ImprovementUtil"+index,"ImprovementUtil"+index,20,20,1,50,95,5,0,times-1);
	SchedulabilityAnalysis::clear_all();
	//DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"ImprovementUtil"+index,"ImprovementUtil"+index,20,20,1,50,50,5,0,times-1);
	//SchedulabilityAnalysis::clear_all();
#endif

	//testMaxHeapMemory();

#if 0
	times = 1000;
	bool deadlineProperty = true;
	//srand (time(NULL));

	for (int fixedUtil = 90; fixedUtil <=90; fixedUtil += 10) {
		//DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2("ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),deadlineProperty,10,25,5,40,5,fixedUtil,fixedUtil,5,0,times-1);
		
		DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),5,40,5,fixedUtil,fixedUtil,5,0,times-1);
		SchedulabilityAnalysis::clear_all();
		DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),5,40,5,fixedUtil,fixedUtil,5,0,times-1);
		SchedulabilityAnalysis::clear_all();
		
	}
#endif

	//DigraphExperiments::testLinearDefectForIBFRBF("Digraphs"+index,result+"TestLinearDefectOfIBFRBF"+index,periodType,0,times-1,1);

	//testRBFIBF("Digraphs"+index,result+"DigraphTimeCalIBF"+index,0,times-1,1);

	//DigraphExperiments::generateDigraphs2("DigraphsTF"+index,50,95,5,times);
	//DigraphExperiments::testDigraphs2("DigraphsTF"+index,result+"DigraphTF"+index, 50,95,5, times);

	//int fixedNodeNum = 2;
	//int fixedNodeNum = 5;
	//int fixedNodeNum = 10;
	//int fixedNodeNum = 20;
	//int fixedNodeNum = 50;
	//DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriod2("DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum),fixedNodeNum,0,times-1);
	//DigraphExperiments::testEffectOfDigraphPeriod("DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum), result+"DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum),0,times-1,3);
	//DigraphExperiments::testEffectOfDigraphPeriod("DigraphPeriod"+index, result+"DigraphPeriod"+index,343,344,1);
#if 0
	times = 1000;
	int util = 10;
	//bool isExact = false;
	int factor = 20;
	srand((unsigned) time(NULL));
	DigraphExperiments::generateDigraphsForTestingEffectOfTaskParameters2("DigraphParameters"+index,40,40,40,5,util,util,10,800,times-1);
	//DigraphExperiments::testEffectOfTaskParameters2("DigraphParameters"+index, result+"DigraphParameters"+index+"Util"+Utility::int_to_string(util)+"factor"+Utility::int_to_string(factor),0.01*factor,40,40,5,util,util,10,933,times-1,1);
#endif

	/*
	DigraphMethodChoice dmc = DigraphRBF;
	DigraphExperiments::generateDigraphsForTestingAcceptanceRatio("DigraphAccRatio"+index,18,60,3,times);
	DigraphExperiments::testAcceptanceRatio(dmc,"DigraphAccRatio"+index,result+"DigraphAccRatio"+index,18,60,3,0,times-1);
	*/

	if (false) {
		//DigraphMethodChoice dmc = DigraphExact;
		//DigraphMethodChoice dmc = DigraphIBFGConstrainedDeadline;
		//DigraphMethodChoice dmc = DigraphIBFLConstrainedDeadline;
		//DigraphMethodChoice dmc = DigraphIBFGLMADDeadline;
		//DigraphMethodChoice dmc = DigraphIBFLLMADDeadline;
		//DigraphMethodChoice dmc = DigraphIBFGArbitraryDeadline;
		//DigraphMethodChoice dmc = DigraphIBFLArbitraryDeadline;
		//DigraphMethodChoice dmc = DigraphRBFG;
		//DigraphMethodChoice dmc = DigraphRBFL;

		bool deadlineProperty = true;
		srand (time(NULL));

#if 0	
		if (argc != 9) {
			cerr << "Error parameter number => " << argc << endl;
			cerr << "Info: *.exe 20 20 50 95 0 999 Util All" << endl;
			cerr << "      *.exe 5 40 60 60 0 999 Num All" << endl;
			exit(EXIT_FAILURE);
		}

		int minTaskNum = atoi(argv[1]);
		int maxTaskNum = atoi(argv[2]);
		int minUtil = atoi(argv[3]);
		int maxUtil = atoi(argv[4]);
		int minRun = atoi(argv[5]);
		int maxRun = atoi(argv[6]); 
		string testName(argv[7]);
		string fileName(argv[8]); 

		if (testName == "Util") {
			DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"Improvement"+testName+index+"file"+fileName,"Improvement"+testName+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"Improvement"+testName+index+"file"+fileName,"Improvement"+testName+index,minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
		}
		else {
			DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"Improvement"+testName+index+"Util"+Utility::int_to_string(minUtil)+"file"+fileName,"Improvement"+testName+index+"Util"+Utility::int_to_string(minUtil),minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
			DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"Improvement"+testName+index+"Util"+Utility::int_to_string(minUtil)+"file"+fileName,"Improvement"+testName+index+"Util"+Utility::int_to_string(minUtil),minTaskNum,maxTaskNum,5,minUtil,maxUtil,5,minRun,maxRun);
			SchedulabilityAnalysis::clear_all();
		}
#endif

#if 1
		if (true) { // test Util
			times = 100;
			DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2("ImprovementUtil"+index,deadlineProperty, 15,20,20,1,50,95,5,0,times-1);
			DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"ImprovementUtil"+index,"ImprovementUtil"+index,20,20,1,50,95,5,0,times-1);
			SchedulabilityAnalysis::clear_all();
			DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"ImprovementUtil"+index,"ImprovementUtil"+index,20,20,1,50,95,5,0,times-1);
			SchedulabilityAnalysis::clear_all();
			//DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2("ImprovementUtil"+index,40,20,20,1,5,5,5,0,0);
			//DigraphExperiments::schedAnalysis(dmc, result+"ImprovementUtil"+index,"ImprovementUtil"+index,20,20,1,90,90,5,650,650);
		} else { // test Num
			int fixedUtil = 60;
			DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2("ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),deadlineProperty,40,40,40,5,fixedUtil,fixedUtil,5,0,times-1);
			DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),10,40,5,fixedUtil,fixedUtil,5,0,times-1);
			SchedulabilityAnalysis::clear_all();
			DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),5,40,5,fixedUtil,fixedUtil,5,0,times-1);
			SchedulabilityAnalysis::clear_all();
		}
#endif	
	}

	//testRBFIBF("TestRBFIBF",result+"TestRBFIBF",0,0,1);

#if 0 // with the fixed vertex number
	int fixedNodeNumbers[5] = {2,5,10,20,50};
	for (int i = 0; i < 5; i ++) {
		int fixedNodeNum = fixedNodeNumbers[i];
		DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriod2("DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum),fixedNodeNum,0,times-1);
		DigraphExperiments::testEffectOfDigraphPeriod("DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum), result+"DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum),0,times-1,3);
	}
#endif

#if 0 // with the fixed vertex number and util
	int fixedNodeNumbers2[2] = {2,20};
	times = 1000;
	for (int i = 0; i < 2; i ++) {
		int fixedNodeNum = fixedNodeNumbers2[i];
		DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriodUtil("DigraphPeriodUtil"+index+"Node"+Utility::int_to_string(fixedNodeNum),fixedNodeNum,5,95,5,0,times-1);
		DigraphExperiments::testEffectOfDigraphPeriodUtil("DigraphPeriodUtil"+index+"Node"+Utility::int_to_string(fixedNodeNum), result+"DigraphPeriodUtil"+index+"Node"+Utility::int_to_string(fixedNodeNum),5,95,5,0,times-1);
	}
#endif

#if 0 // with the fixed utilization and factor
	int utils[3] = {30,50,70}; 
	int factors[3] = {20,50,80};

	bool isExact = false;
	for (int i=0; i<3; i++) {
		int util = utils[i];
		DigraphExperiments::generateDigraphsForTestingEffectOfTaskParameters2("DigraphParameters"+index,40,5,40,5,util,util,10,0,times-1);
		for (int j=0; j<3; j++) {
			int factor = factors[j];
			DigraphExperiments::testEffectOfTaskParameters2("DigraphParameters"+index, result+"DigraphParameters"+index+"Util"+Utility::int_to_string(util)+"factor"+Utility::int_to_string(factor),0.01*factor,5,40,5,util,util,10,0,times-1,1);
		}
	}
#endif

#if 0 // with the fixed task and factor
	times = 1000;

	int fixedTaskNumber[2] = {10,30};
	int factors[3] = {20,50,80};

	for (int i=0; i<2; i++) {
		int taskNumber = fixedTaskNumber[i];
		DigraphExperiments::generateDigraphsForTestingEffectOfTaskParameters2("DigraphParametersTaskNumber"+index,40,taskNumber,taskNumber,5,10,90,10,0,times-1);
		for (int j=0; j<3; j++) {
			int factor = factors[j];
			DigraphExperiments::testEffectOfTaskParameters2Util("DigraphParametersTaskNumber"+index, result+"DigraphParameters"+index+"TaskNum"+Utility::int_to_string(taskNumber)+"factor"+Utility::int_to_string(factor),0.01*factor,taskNumber,taskNumber,5,10,90,10,0,times-1,1);
		}
	}
#endif

	return 0;
}