#include "FSMExperiments.h"

string FSMExperiments::funcNames[28] = {
	"EXACTSO", // Exact schedulability analysis on FSMs with static offsets
	"EXACTAO1", // Exact schedulability analysis on (precise) digraphs transformed with action instances
	"EXACTAO2", // Exact schedulability analysis on (simple) digraphs transfromed with actions
	"IBFSO", // IBF schedulability analysis on FSMs with static offsets
	"RBFSO", // RBF schedulability analysis on FSMs with static offsets 
	"IBFSOP", // IBF schedulability analysis on FSMs with static offsets using the linear periodicity
	"RBFSOP", // RBF schedulability analysis on FSMs with static offsets using the linear periodicity
	"IBFAOFSM", // IBF schedulability analysis with arbitrary offsets on FSMs
	"RBFAOFSM", // RBF schedulability analysis with arbitrary offsets on FSMs
	"IBFAOFSMDP", // IBF schedulability analysis with arbitrary offsets on FSMs using dynamic programming techniques
	"RBFAOFSMDP", // RBF schedulability analysis with arbitrary offsets on FSMs using dynamic programming techniques
	"IBFAO1G", // IBF schedulability analysis on (precise) digraphs using Guan's algorithm
	"RBFAO1G", // RBF schedulability analysis on (precise) digraphs using Guan's algorithm
	"IBFAO1DP", // IBF schedulability analysis on (precise) digraphs using dynamic programming techniques
	"RBFAO1DP", // RBF schedulability analysis on (precise) digraphs using dynamic programming techniques
	"IBFAO1L", // IBF schedulability analysis on (precise) digraphs using the linear periodicity
	"RBFAO1L", // RBF schedulability analysis on (precise) digraphs using the linear periodicity
	"IBFAO1LSCC", // IBF schedulability analysis on (strongly-connected precise) digraphs using the linear periodicity
	"RBFAO1LSCC", // RBF schedulability analysis on (strongly-connected precise) digraphs using the linear periodicity
	"IBFAO2G", // IBF schedulability analysis on (simple) digraphs using Guan's algorithm
	"RBFAO2G", // RBF schedulability analysis on (simple) digraphs using Guan's algorithm
	"IBFAO2DP", // IBF schedulability analysis on (simple) digraphs using dynamic programming techniques
	"RBFAO2DP", // RBF schedulability analysis on (simple) digraphs using dynamic programming techniques
	"IBFAO2L", // IBF schedulability analysis on (simple) digraphs using the linear periodicity
	"RBFAO2L", // RBF schedulability analysis on (simple) digraphs using the linear periodicity
	"LUBIBF", // LUBIBF schedulability analysis on FSMs with static offsets
	"LUBRBF", // LUBRBF schedulability analysis on FSMs with static offsets
	"LUBCSUM" // LUBCSUM schedulability analysis on FSMs with static offsets
};

void FSMExperiments::generateRandomSystemsForExactAnalysis(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int minRun, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;

	int maxState = 15;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int run=minRun; run<=maxRun; run++) {
					try {
						Stateflow** sfs;
						while(true) {
							sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis7(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);
							break;
							//SchedulabilityAnalysis::generate_critical_action_pair(sfs,num);
							//if (SchedulabilityAnalysis::generate_request_function(sfs,num,false)) break;

							for (int i=0; i<num; i++)
								delete sfs[i];
							delete[] sfs;
						}

						string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						cout<<name<<endl;

						FileWriter::DotFileWriter(sfs, num, p);

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;

					}
					catch(bad_alloc) // catch bad_alloc exception
					{
						cout<<"Exception"<<endl;
					}
				}
			}
		}
	}
}

void FSMExperiments::generateRandomSystemsForApproximateAnalysis(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;

	int maxState = 15;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int run=0; run<maxRun; run++) {
					try {
						//Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);
						Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis4(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);

						string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						cout<<name<<endl;

						FileWriter::DotFileWriter(sfs, num, p);

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;

					}
					catch(bad_alloc) // catch bad_alloc exception
					{
						cout<<"Exception"<<endl;
					}
				}
			}
		}
	}
}

void FSMExperiments::generateRandomSystemsForApproximateAnalysisWithPeriodicity(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;

	int maxState = 15;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int run=0; run<maxRun; run++) {
					try {
						//Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);
						Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis_with_periodicity(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);

						string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						cout<<name<<endl;

						FileWriter::DotFileWriter(sfs, num, p);

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;

					}
					catch(bad_alloc) // catch bad_alloc exception
					{
						cout<<"Exception"<<endl;
					}
				}
			}
		}
	}
}

void FSMExperiments::generateRandomSystemsForApproximateAnalysisWithPeriodicity2(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int minRun, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;

	int maxState = 25;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int run=minRun; run<=maxRun; run++) {
					try {
						//Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);
						Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis_with_periodicity(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);

						string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
						const char *p = name.c_str();

						cout<<name<<endl;

						FileWriter::DotFileWriter(sfs, num, p);

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;

					}
					catch(bad_alloc) // catch bad_alloc exception
					{
						cout<<"Exception"<<endl;
					}
				}
			}
		}
	}
}

void FSMExperiments::SchedAnalysis(MethodChoice mc, string file, string directory,int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int startRun, int endRun) {
	// write a file
	string result = file+funcNames[mc];
	const char *p = result.c_str();

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int choice = 0;
	int maxState = 15;
	bool output = false;

	SchedulabilityAnalysis::set_zero();
	map<int,map<int,vector<double>>> allRecTimes;
	Timer timer;

	for (int num=minNum; num<=maxNum; num+=stepNum) {

		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {

			for (int scc = minScc; scc <= maxScc; scc += stepScc) {
				int nException = 0;
				int nUnSuccess = 0;

				vector<double> recTimes;
				double sumTUtil2 = 0;

				for (int run=startRun; run<=endRun; run++) {
					Stateflow** sfs;

					int numStateflows;
					string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
					const char* file = name.c_str();
					FileReader::DotFileReader(sfs, numStateflows,1,file);

					try {
						timer.start();
						bool schedulable = doSchedAnalysis(mc,sfs,num,fout);
						timer.end();
						recTimes.push_back(timer.getTime());

						cout<<file<<"=>"<< FSMExperiments::funcNames[mc] << "\t" << schedulable << endl;
						fout<<file<<"=>"<< FSMExperiments::funcNames[mc] << "\t" << schedulable << endl;

						double tUtil = 0;
						double tUtil2 = 0;
						for (int i=0; i<num; i++) {
							Stateflow* sf = sfs[i];
							tUtil+=sf->lfac;

							if (sf->isIrred) SchedulabilityAnalysis::nIrreducibleStateflows++;
							/*
							if (sf->simple_digraph != NULL)
								tUtil2+=sf->simple_digraph->linear_factor;
								*/
						}

						sumTUtil2 += tUtil2;

						cout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil <<", Actual utilization = "<<tUtil2<<endl;
						fout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil <<", Actual utilization = "<<tUtil2<<endl;

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
					} 
					catch(bad_alloc& ba)
					{
						nException++;

						for (int i=0; i<num; i++)
							delete sfs[i];
						delete[] sfs;
					}
				}

				allRecTimes[num][util] = recTimes;

				fout<<"nException="<<nException<<endl;

				SchedulabilityAnalysis::nStateflows = num;
				SchedulabilityAnalysis::totalUtilization = util;
				SchedulabilityAnalysis::totalUtilization2 = sumTUtil2/(endRun-startRun+1);
				// reset
				SchedulabilityAnalysis::reset();

				fout<<"Chocice="<<choice<<endl;

				SchedulabilityAnalysis::output_vectors(fout);
				fout << "=============================================" << endl;
				FSMExperiments::output(allRecTimes,funcNames[mc],fout);
			}
		}
	}
	fout<<"Chocice="<<choice<<endl;

	SchedulabilityAnalysis::output_vectors(fout);
	fout << "=============================================" << endl;
	FSMExperiments::output(allRecTimes,funcNames[mc],fout);

	fout.close();
}

void FSMExperiments::generateRandomSystemsForWeightedSchedAnalysis(string directory, int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int minScc, int maxScc, int stepScc, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;

	int maxState = 15;
	for (int num=minNum; num<=maxNum; num+=stepNum) {
		//for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
		for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
			for (int run=0; run<maxRun; run++) {
				try {
					int util = rand()%(maxUtil-minUtil)+minUtil; // util \in [minUtil,maxUtil]
					Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis_with_periodicity(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);

					string name = directory+"\\WSNum"+Utility::int_to_string(num)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
					const char *p = name.c_str();

					cout<<name<<endl;

					FileWriter::DotFileWriter(sfs, num, p);

					for (int i=0; i<num; i++)
						delete sfs[i];
					delete[] sfs;

				}
				catch(bad_alloc) // catch bad_alloc exception
				{
					cout<<"Exception"<<endl;
				}
			}
		}
		//}
	}
}

void FSMExperiments::WeightedSchedulabilityAnalysis(MethodChoice mc,string directory, string file, int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int minScc, int maxScc, int stepScc, int maxRun) {
	// #define PERIODICITY_PROPERTY
	// write a file
	string result = file+funcNames[mc];
	const char *p = result.c_str();

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int choice = 2;
	int maxState = 15;
	bool output = false;

	vector<double> vec_sum_weightedSched;
	
	SchedulabilityAnalysis::set_zero();

	for (int num=minNum; num<=maxNum; num+=stepNum) {

		//for (int util = minUtil; util<=maxUtil; util+=stepUtil) {

		for (int scc = minScc; scc <= maxScc; scc += stepScc) {
			int nException = 0;

			double sum_weightedSched = 0;

			double totalUtil = 0;

			for (int run=0; run<maxRun; run++) {
				Stateflow** sfs;

				int numStateflows;
				string name = directory+"\\WSNum"+Utility::int_to_string(num)+"Scc"+Utility::int_to_string(scc)+"Run"+Utility::int_to_string(run)+".dot";
				const char* file = name.c_str();
				FileReader::DotFileReader(sfs, numStateflows,1,file);

				try {
					bool schedulable = doSchedAnalysis(mc,sfs,num,fout);

					cout<<file<<"=>"<< FSMExperiments::funcNames[mc] << "\t" << schedulable << endl;
					fout<<file<<"=>"<< FSMExperiments::funcNames[mc] << "\t" << schedulable << endl;

					if (mc != LUBIBF && mc != LUBRBF && mc != LUBCSUM) {
						// we need to calculate the util for the FSMs
						for (int i=0; i<num; i++) {
							Stateflow* sf = sfs[i];
							sf->generate_rbfs();
							sf->calculate_linear_factor();
						}
					}

					double tUtil = 0;
					for (int i=0; i<num; i++) {
						Stateflow* sf = sfs[i];
						tUtil+=sf->lfac;
					}

					totalUtil += tUtil;

					sum_weightedSched += tUtil * schedulable;

					for (int i=0; i<num; i++)
						delete sfs[i];
					delete[] sfs;
				} 
				catch(bad_alloc& ba)
				{
					nException++;

					for (int i=0; i<num; i++)
						delete sfs[i];
					delete[] sfs;
				}
			}
			fout<<"nException="<<nException<<endl;
			SchedulabilityAnalysis::nStateflows = num;

			vec_sum_weightedSched.push_back(sum_weightedSched/totalUtil);

			//SchedulabilityAnalysis::totalUtilization = util;
			// reset
			SchedulabilityAnalysis::reset();
			SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSched"+FSMExperiments::funcNames[mc],vec_sum_weightedSched);
		}
		//}
	}
	fout<<"Chocice="<<choice<<endl;

	SchedulabilityAnalysis::output_vectors(fout);
	SchedulabilityAnalysis::output_one_vector(fout,"sum_weightedSched"+FSMExperiments::funcNames[mc],vec_sum_weightedSched);

	fout.close();
}

void FSMExperiments::generateRandomSystemsForScalabilityTest(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int minState, int maxState, int stepState, int minRun,int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	RandomGenerator::STATEFLOW_SCALE = 10000;

	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			for (int scc = minScc; scc<=maxScc; scc+=stepScc) {
				for (int state = minState; state <= maxState; state += stepState) { 
					for (int run=minRun; run<=maxRun; run++) {
						try {
							//Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis(num, maxState, 1.0*util/100, period_choice,1.0*scc/100);
							Stateflow** sfs = RandomGenerator::generate_stateflow_system_for_approximate_analysis3(num, state, 1.0*util/100, period_choice,1.0*scc/100);

							string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"State"+Utility::int_to_string(state)+"Run"+Utility::int_to_string(run)+".dot";
							const char *p = name.c_str();

							cout<<name<<endl;

							FileWriter::DotFileWriter(sfs, num, p);

							for (int i=0; i<num; i++)
								delete sfs[i];
							delete[] sfs;

						}
						catch(bad_alloc) // catch bad_alloc exception
						{
							cout<<"Exception"<<endl;
						}
					}
				}
			}
		}
	}
}

void FSMExperiments::SchedAnalysisForScalabilityTest(MethodChoice mc, bool myProperty, string file, string directory, int choice, int num, int util, int scc, int minState, int maxState, int stepState, int startRun, int endRun) {
	// write a file
	string result = file+funcNames[mc];
	const char *p = result.c_str();

	ofstream fout(result, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<result<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	bool output = false;

	SchedulabilityAnalysis::set_zero();

	for (int state = minState; state <= maxState; state += stepState) {
		int nException = 0;
		int nUnSuccess = 0;
		int nTimeout = 0;

		for (int run=startRun; run<=endRun; run++) {
			Stateflow** sfs;

			int numStateflows;
			string name = directory+"\\Num"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Scc"+Utility::int_to_string(scc)+"State"+Utility::int_to_string(state)+"Run"+Utility::int_to_string(run)+".dot";
			const char* file = name.c_str();
			FileReader::DotFileReader(sfs, numStateflows,1,file);

			try {
				bool schedulable = doSchedAnalysis(mc,sfs,num,fout);
				
				cout<<file<<"=>"<< FSMExperiments::funcNames[mc] << "\t" << schedulable << endl;
				fout<<file<<"=>"<< FSMExperiments::funcNames[mc] << "\t" << schedulable << endl;

				double tUtil = 0;
				for (int i=0; i<num; i++) {
					Stateflow* sf = sfs[i];
					tUtil+=sf->lfac;
				}
				cout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil<<endl;
				fout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil<<endl;

				for (int i=0; i<num; i++)
					delete sfs[i];
				delete[] sfs;
			} 
			catch(bad_alloc& ba)
			{
				nException++;

				for (int i=0; i<num; i++)
					delete sfs[i];
				delete[] sfs;
			}
		}
		fout<<"nException="<<nException<<endl;

		SchedulabilityAnalysis::nStateflows = num;
		SchedulabilityAnalysis::totalUtilization = util;

		// reset
		SchedulabilityAnalysis::reset();

		fout<<"Chocice="<<choice<<endl;

		SchedulabilityAnalysis::output_vectors(fout);
	}

	fout<<"Chocice="<<choice<<endl;

	SchedulabilityAnalysis::output_vectors(fout);

	fout.close();
}

bool FSMExperiments::doSchedAnalysis(MethodChoice mc, Stateflow** sfs, int num, ostream& out) {
	bool schedulable;
	switch (mc)
	{
	case EXACTSO:
		schedulable = SchedulabilityAnalysis::exact_sched_analysis(sfs,num,0,false,out);
		if (schedulable) SchedulabilityAnalysis::nExactStaticOffset++;
		break;
	case EXACTAO1:
		schedulable = SchedulabilityAnalysis::exact_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,0);
		//schedulable = SchedulabilityAnalysis::exact_analysis_arbitrary_offset_by_precise_digraphs_timeout(sfs,num,0,10);
		if (schedulable) SchedulabilityAnalysis::nExactArbitraryOffsetByPreciseDigraph++;
		break;
	case EXACTAO2:
		schedulable = SchedulabilityAnalysis::exact_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,0);
		if (schedulable) SchedulabilityAnalysis::nExactArbitraryOffsetBySimpleDigraph++;
		break;
	case IBFSO:
		schedulable = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nIBFStaticOffset++;
		break;
	case RBFSO:
		schedulable = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nRBFStaticOffset++;
		break;
	case IBFSOP:
		schedulable = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nIBFStaticOffset++;
		break;
	case RBFSOP:
		schedulable = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nRBFStaticOffset++;
		break;
	case RBFAOFSM:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(sfs,num,0);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffset++;
		break;
	case IBFAOFSM:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(sfs,num,0);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffset++;
		break;
	case RBFAOFSMDP:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_reachability_digraphs_DP(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffset++;
		break;
	case IBFAOFSMDP:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_reachability_digraphs_DP(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffset++;
		break;
	case IBFAO1G:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph++;
		break;
	case RBFAO1G:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph++;
		break;
	case IBFAO1DP:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs_DP(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph++;
		break;
	case RBFAO1DP:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs_DP(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph++;
		break;
	case IBFAO1L:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs_DP(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraphWithLinearization++;
		break;
	case RBFAO1L:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs_DP(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraphWithLinearization++;
		break;
	case IBFAO1LSCC:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_strongly_connected_precise_digraphs_DP(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraphWithLinearization++;
		break;
	case RBFAO1LSCC:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_strongly_connected_precise_digraphs_DP(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraphWithLinearization++;
		break;
	case IBFAO2G:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph++;
		break;
	case RBFAO2G:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph++;
		break;
	case IBFAO2DP:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs_DP(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph++;
		break;
	case RBFAO2DP:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs_DP(sfs,num,0,false);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph++;
		break;
	case IBFAO2L:
		schedulable = SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraphWithLinearization++;
		break;
	case RBFAO2L:
		schedulable = SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs_DP(sfs,num,0,true);
		if (schedulable) SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraphWithLinearization++;
		break;
	case LUBIBF:
		schedulable = SchedulabilityAnalysis::lu_ibf_sched_analysis(sfs,num,2);
		if (schedulable) SchedulabilityAnalysis::nLinearUpperIBF++;
		break;
	case LUBRBF:
		schedulable = SchedulabilityAnalysis::lu_rbf_sched_analysis(sfs,num,2);
		if (schedulable) SchedulabilityAnalysis::nLinearUpperRBF++;
		break;
	case LUBCSUM:
		schedulable = SchedulabilityAnalysis::lu_csum_sched_analysis(sfs,num,0);
		if (schedulable) SchedulabilityAnalysis::nLinearUpperCSUM++;
		break;
	default:
		cerr << "Error the method choice: " << mc << endl;
		exit(EXIT_FAILURE);
		break;
	}
	return schedulable;
}

void FSMExperiments::output(map<int,map<int,vector<double>>> allRecTimes, string name, ostream& out) {
	for (map<int,map<int,vector<double>>>::iterator iter0 = allRecTimes.begin(); iter0 != allRecTimes.end(); iter0++) {
		for (map<int,vector<double>>::iterator iter1 = iter0->second.begin(); iter1 != iter0->second.end(); iter1++) {
			out << name << "[" << iter0->first <<"][" << iter1->first << "]=[";

			for (vector<double>::iterator iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
				out << *iter2;
				if (iter2 != (--iter1->second.end()))
					out << ",";
			}

			out << "];" << endl;
			out << endl;
		}

	}
}