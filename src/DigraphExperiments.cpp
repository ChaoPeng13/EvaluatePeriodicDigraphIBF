#include "DigraphExperiments.h"

using namespace std;

string DigraphExperiments::digraphFuncNames[11] = {
	"DigraphExact",
	"DigraphIBFGConstrainedDeadline",
	"DigraphIBFLConstrainedDeadline",
	"DigraphIBFGLMADDeadline",
	"DigraphIBFLLMADDeadline",
	"DigraphIBFGArbitraryDeadline",
	"DigraphIBFLArbitraryDeadline",
	"DigraphRTIBFGConstrainedDeadline",
	"DigraphRTIBFLConstrainedDeadline",
	"DigraphRBFG",
	"DigraphRBFL"
};

bool DigraphExperiments::doDigraphSchedAnalysis(DigraphMethodChoice dmc, Digraph** digraphs, int num) {
	bool schedulable;
	switch (dmc)
	{
	case DigraphExact:
		schedulable = SchedulabilityAnalysis::digraph_exact_analysis(digraphs,num,-1);
		if (schedulable) SchedulabilityAnalysis::nDigraphExact++;
		break;
	case DigraphIBFGConstrainedDeadline:
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(digraphs,num,-1,false,false);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphIBFLConstrainedDeadline:
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(digraphs,num,-1,true,false);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphIBFGLMADDeadline:
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(digraphs,num,-1,false,true);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphIBFLLMADDeadline:
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(digraphs,num,-1,true,true);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphIBFGArbitraryDeadline:
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs,num,-1,false);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphIBFLArbitraryDeadline:
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs,num,-1,true);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphRTIBFGConstrainedDeadline:
		schedulable = SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_constrained_deadline(digraphs,num,-1,false);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphRTIBFLConstrainedDeadline:
		schedulable = SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_constrained_deadline(digraphs,num,-1,true);
		if (schedulable) SchedulabilityAnalysis::nDigraphIBF++;
		break;
	case DigraphRBFG:
		schedulable = SchedulabilityAnalysis::digraph_rbf_analysis_dynamic_programming(digraphs,num,-1,false,false);
		if (schedulable) SchedulabilityAnalysis::nDigraphRBF++;
		break;
	case DigraphRBFL:
		schedulable = SchedulabilityAnalysis::digraph_rbf_analysis_dynamic_programming(digraphs,num,-1,true,false);
		if (schedulable) SchedulabilityAnalysis::nDigraphRBF++;
		break;
	default:
		cerr << "Error the method choice: " << dmc << endl;
		exit(EXIT_FAILURE);
		break;
	}
	return schedulable;
}

void DigraphExperiments::generateDigraphsForTestingAcceptanceRatio(string directory, int minUtil, int maxUtil, int stepUtil, int totalRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		std::cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	for (int util = minUtil; util <= maxUtil; util += stepUtil) {
		for (int run = 0; run <totalRun; run++) {
			int numTask;
			Digraph** digraphs = RandomGenerator::generate_digraph_system_guan(numTask,0.01*util);

			double realUtil = 0;
			for (int i=0; i<numTask; i++) {
				Digraph* digraph = digraphs[i];
				realUtil += digraph->linear_factor;
			}
			cout << "ExpectedUtil=" << 0.01*util << "\tRealUtil=" << realUtil << endl;
			
			string name = directory+"\\digraphUtil"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run);
			const char *p = name.c_str();
			cout << name <<endl;
			FileWriter::DotFileWriter(digraphs, numTask, p);

			for(int i=0; i<numTask; i++)
				delete digraphs[i];
			delete[] digraphs;
		}
	}
}

void DigraphExperiments::testAcceptanceRatio(DigraphMethodChoice dmc, string directory, string file, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun) {
	file = file + digraphFuncNames[dmc];

	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	bool output = false;

	SchedulabilityAnalysis::set_zero();

	for (int util = minUtil; util <= maxUtil; util += stepUtil) {
		int nException = 0;
		for (int run=startRun; run <= endRun; run++) {
			Digraph** digraphs;

			int numDigraphs;
			string name = directory+"\\digraphUtil"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run);
			const char *p = name.c_str();
			FileReader::DotFileReader(digraphs,numDigraphs,1,p);

			try {
				bool schedulable = doDigraphSchedAnalysis(dmc,digraphs,numDigraphs);
				cout<<name<<"=>"<< DigraphExperiments::digraphFuncNames[dmc] << "\t" << schedulable << endl;
				fout<<name<<"=>"<< DigraphExperiments::digraphFuncNames[dmc] << "\t" << schedulable << endl;

				for (int i=0; i<numDigraphs; i++)
					delete digraphs[i];
				delete[] digraphs;
			} 
			catch(bad_alloc& ba)
			{
				nException++;

				for (int i=0; i<numDigraphs; i++)
					delete digraphs[i];
				delete[] digraphs;
			}
		}
		fout<<"nException="<<nException<<endl;

		SchedulabilityAnalysis::totalUtilization = util;

		// reset
		SchedulabilityAnalysis::reset();
		SchedulabilityAnalysis::output_vectors(fout);
	}

	SchedulabilityAnalysis::output_vectors(fout);
	fout.close();
}

void DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation(string directory) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		std::cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int numTask = 20;
	int run = 10;
	for (int util = 10; util <= 90; util += 10) {
		int totalNum = 0;
		while (totalNum < run) {
			Digraph** digraphs = RandomGenerator::generate_digraph_system3(numTask,0.01*util);

			double realUtil = 0;
			for (int i=0; i<numTask; i++) {
				Digraph* digraph = digraphs[i];
				realUtil += digraph->linear_factor;
			}
			cout << "ExpectedUtil=" << 0.01*util << "\tRealUtil=" << realUtil << endl;
			if (abs(0.01*util-realUtil) <= 0.01) {
				cout << "found\t" << totalNum << endl;

				string name = directory+"\\digraphUtil"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(totalNum);
				const char *p = name.c_str();
				cout << name <<endl;
				FileWriter::DotFileWriter(digraphs, numTask, p);

				totalNum++;
			}

			delete[] digraphs;
		}
	}
}

void DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation(string directory,  int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int totalRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	for (int numTask = minNum; numTask <= maxNum; numTask += stepNum) {
		for (int util = minUtil; util <= maxUtil; util += stepUtil) {
			int index = 0;
			while (index < totalRun) {
				Digraph** digraphs = RandomGenerator::generate_mixed_digraph_system(numTask,0.01*util);

				double realUtil = 0;
				for (int i=0; i<numTask; i++) {
					Digraph* digraph = digraphs[i];
					realUtil += digraph->linear_factor;
				}

				cout << "ExpectedUtil=" << 0.01*util << "\tRealUtil=" << realUtil << endl;

				if (abs(0.01*util-realUtil) <= 0.01) {
					string name = directory+"\\digraphNum"+Utility::int_to_string(numTask)+"Util"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(index)+".dot";
					const char *p = name.c_str();
					cout << name <<endl;
					FileWriter::DotFileWriter(digraphs, numTask, p);
					index++;
				}

				//digraphs[0]->write_graphviz(cout);

				for (int i=0; i<numTask; i++) {
					//cout << "Relese digraph " << i <<endl;
					delete digraphs[i];
				}
				delete[] digraphs;
			}
		}
	}
}

void DigraphExperiments::testImprovementOnIBFCalculationInDigraphSystem(string directory, string file, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int totalRun) {
	const char* p = file.c_str();

	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	SchedulabilityAnalysis::set_zero();

	vector<double> vec_imp;
	vector<double> vec_caltime_with_periodicity;
	vector<double> vec_caltime_without_periodicity;
	vector<int> vec_success_num; // r < tf;

	Timer timer;
	for (int numTask = minNum; numTask <= maxNum; numTask += stepNum) {
		for (int util = minUtil; util <= maxUtil; util += stepUtil) {
			double imp = 0;
			double caltime_with_periodicity = 0;
			double caltime_without_periodicity = 0;
			int success_num = 0;

			for (int run=0; run<totalRun; run++) {
				//string input = directory+"\\digraphUtil"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run);
				string input = directory+"\\digraphNum"+Utility::int_to_string(numTask)+"Util"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run)+".dot";
				const char *p = input.c_str();
				cout << input <<endl;
				fout << "Read file " << input << endl;
				Digraph** digraphs;
				int num;
				FileReader::DotFileReader(digraphs,num, 1, p);

				SchedulabilityAnalysis::prepare_all_digraphs2(digraphs,num);

				for (int i=0; i<num; i++) {
					Digraph* digraph = digraphs[i];
					if (i==num-1) {
						cout << "Improvement=" << (digraph->tf1-digraph->tf2)/digraph->tf1 << endl;
						fout << "Improvement=" << (digraph->tf1-digraph->tf2)/digraph->tf1 << endl;
						imp += (digraph->tf1-digraph->tf2)/digraph->tf1;
					}
				}

				// do the schedulability analysis
				SchedulabilityAnalysis::generate_critical_vertices(digraphs,num);

				/*
				// set tf, denoted as the maximum time length for calculating ibf
				int tf_max = INT_MIN;
				for (int i=0; i<num; i++) {
					Digraph* digraph = digraphs[i];
					for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
						tf_max = max(tf_max,(*iter)->deadline);
					}
				}

				cout << "tf_max=" << tf_max << endl;
				*/
			
				double sumTime0 = 0;
				double sumTime1 = 0;

				/*
				double sumCRBF=0, sumCIBF=0, sumUtil=0;
				for (int i=0; i<num-1; i++) {
					Digraph* digraph = digraphs[i];
					cout << "tf=" << digraph->tf << "\tlfac=" << digraph->linear_factor << endl;
					cout << "tf0=" << digraph->tf0 << "\ttf1=" << digraph->tf1 << "\ttf2=" << digraph->tf2 << endl;
					sumCRBF += digraph->c_rbf;
					sumCIBF += digraph->c_ibf;
					sumUtil += digraph->linear_factor;
				}

				sumCRBF += digraphs[num-1]->c_dbf;
				sumCIBF += digraphs[num-1]->c_dbf;
				sumUtil += digraphs[num-1]->linear_factor;
				cout << "sumCRBF=" << sumCRBF << "\tsumCIBF=" << sumCIBF << "\tsumUtil=" << sumUtil << endl;
				*/
				cout << "tf0=" << digraphs[num-1]->tf0 << "\ttf1=" << digraphs[num-1]->tf1 << "\ttf2=" << digraphs[num-1]->tf2 << endl;

				for (int i=0; i<num-1; i++) {
					Digraph* digraph = digraphs[i];
					digraph->tf = digraphs[num-1]->tf2;
					//cout << "tf=" << digraph->tf << "\tlfac=" << digraph->linear_factor << endl;
					//cout << "tf0=" << digraph->tf0 << "\ttf1=" << digraph->tf1 << "\ttf2=" << digraph->tf2 << endl;
					if (digraph->pGCD == 1) {
						double time0=0,time1=0; 
						cout << "Processing digraph " << i << " in " << input << endl;
						cout << "Calculating ibf with periodicity..." << endl;
						map<int,double> ibfs;
						timer.start();
						bool success = digraph->calculate_ibf_with_periodicity(false,ibfs);
						timer.end();
						time0 = timer.getTime();

						SchedulabilityAnalysis::tDigraphCalLinearPeriod += digraph->gran_digraph->tCalLinearPeriod;
						SchedulabilityAnalysis::tDigraphCalLinearDefect += digraph->gran_digraph->tCalLinearDefect;

						cout << "Calculating ibf without periodicity..." << endl;
						cout << digraph->tf << endl;
						timer.start();
						digraph->calculate_ibf_without_periodicity(digraph->tf, digraph->ibf_map2);
						timer.end();
						time1 = timer.getTime();
					

						cout << "success=" << success << "\ttime0=" << time0 << "\ttime1=" << time1 << endl;
						fout << "success=" << success << "\ttime0=" << time0 << "\ttime1=" << time1 << endl;

						// now we test these generated ibfs
						if (success) {
							success_num ++;
							sumTime0 += time0;
							sumTime1 += time1;
							if (true) {
								for (int t=0; t<=digraph->tf; t++) {
									if (abs(digraph->ibf_map3[t]-digraph->ibf_map2[t]) > Utility::EPSILON ) {
										cerr << "Error ibf on t=" << t <<"\t ibf with periodicity=" << digraph->ibf_map3[t] 
											 <<"\tibf without periodicity=" << digraph->ibf_map2[t] << "\n"
											 << "Periodicity calculating: ibf(t-lper)=" << digraph->ibf_map2[t-digraph->gran_digraph->lper]
											 << "\twPeriodicity calculating: ibf(t-lper)=" << digraph->ibf_map3[t-digraph->gran_digraph->lper] << "\n"
											 << "\tLinearFactor=" << digraph->linear_factor 
											 << "\tLinearPeriod=" << digraph->gran_digraph->lper
											 << "\tUbLinearDefect=" << digraph->gran_digraph->ubldef << endl;
										digraph->write_graphviz(cout);
										exit(EXIT_FAILURE);
									}
								}
							}
						}
					}
					else
						digraph->calculate_ibf_without_periodicity(digraphs[num-1]->tf2, digraph->ibf_map2);
				}
				cout << "sumTime0=" << sumTime0 << "\tsumTime1=" << sumTime1 << endl;
				fout << "sumTime0=" << sumTime0 << "\tsumTime1=" << sumTime1 << endl;

				caltime_with_periodicity += sumTime0;
				caltime_without_periodicity += sumTime1;

				for (int i=0; i<num; i++)
					delete digraphs[i];
				delete[] digraphs;
			}

			vec_imp.push_back(imp);
			vec_caltime_with_periodicity.push_back(caltime_with_periodicity);
			vec_caltime_without_periodicity.push_back(caltime_without_periodicity);
			vec_success_num.push_back(success_num);

			SchedulabilityAnalysis::reset();
		}
	}

	SchedulabilityAnalysis::output_vectors(fout);
	SchedulabilityAnalysis::output_one_vector(fout,"Improvement",vec_imp);
	SchedulabilityAnalysis::output_one_vector(fout,"tCalTimeWithPeriodicity",vec_caltime_with_periodicity);
	SchedulabilityAnalysis::output_one_vector(fout,"tCalTimeWithoutPeriodicity",vec_caltime_without_periodicity);
	SchedulabilityAnalysis::output_one_vector(fout,"SuccessNumber", vec_success_num);

	fout.close();
}

void DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2(string directory, bool deadlineProperty, int maxNode, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	for (int numTask = minNum; numTask <= maxNum; numTask += stepNum) {
		for (int util = minUtil; util <= maxUtil; util += stepUtil) {
			for (int run = startRun; run<=endRun; run++) {
				Digraph** digraphs;
				try {
					digraphs = RandomGenerator::generate_digraph_system(numTask,maxNode,0.01*util,deadlineProperty);

					double realUtil = 0;
					for (int i=0; i<numTask; i++) {
						Digraph* digraph = digraphs[i];
						realUtil += digraph->linear_factor;
					}

					cout << "ExpectedUtil=" << 0.01*util << "\tRealUtil=" << realUtil << endl;

					string name = directory+"\\digraphNum"+Utility::int_to_string(numTask)+"Util"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run)+".dot";
					const char *p = name.c_str();
					cout << name <<endl;
					FileWriter::DotFileWriter(digraphs, numTask, p);
					
					//digraphs[0]->write_graphviz(cout);

					for (int i=0; i<numTask; i++) {
						//cout << "Relese digraph " << i <<endl;
						delete digraphs[i];
					}
					delete[] digraphs;
				}
				catch(bad_alloc& ba)
				{
					for (int i=0; i<numTask; i++)
						delete digraphs[i];
					delete[] digraphs;
				}
			}
		}
	}
}

void DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2(string directory, bool deadlineProperty, int minNode, int maxNode, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	for (int numTask = minNum; numTask <= maxNum; numTask += stepNum) {
		for (int util = minUtil; util <= maxUtil; util += stepUtil) {
			for (int run = startRun; run<=endRun; run++) {
				Digraph** digraphs;
				try {
					digraphs = RandomGenerator::generate_digraph_system(numTask,minNode,maxNode,0.01*util,deadlineProperty);

					double realUtil = 0;
					for (int i=0; i<numTask; i++) {
						Digraph* digraph = digraphs[i];
						realUtil += digraph->linear_factor;
					}

					cout << "ExpectedUtil=" << 0.01*util << "\tRealUtil=" << realUtil << endl;

					string name = directory+"\\digraphNum"+Utility::int_to_string(numTask)+"Util"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run)+".dot";
					const char *p = name.c_str();
					cout << name <<endl;
					FileWriter::DotFileWriter(digraphs, numTask, p);

					//digraphs[0]->write_graphviz(cout);

					for (int i=0; i<numTask; i++) {
						delete digraphs[i];
					}
					delete[] digraphs;
				}
				catch(bad_alloc& ba)
				{
					cerr << "Bad Allocation!" << endl;
					//exit(EXIT_FAILURE);
					for (int i=0; i<numTask; i++) {
						delete digraphs[i];
					}
					delete[] digraphs;
				}
			}
		}
	}
}

void DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2WithoutPointer(string directory, bool deadlineProperty, int minNode, int maxNode, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	for (int numTask = minNum; numTask <= maxNum; numTask += stepNum) {
		for (int util = minUtil; util <= maxUtil; util += stepUtil) {
			for (int run = startRun; run<=endRun; run++) {
				Digraph** digraphs;
				try {
					digraphs = RandomGenerator::generate_digraph_system_without_pointer(numTask,minNode,maxNode,0.01*util,deadlineProperty);

					double realUtil = 0;
					for (int i=0; i<numTask; i++) {
						Digraph* digraph = digraphs[i];
						realUtil += digraph->linear_factor;
					}

					cout << "ExpectedUtil=" << 0.01*util << "\tRealUtil=" << realUtil << endl;

					string name = directory+"\\digraphNum"+Utility::int_to_string(numTask)+"Util"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run)+".dot";
					const char *p = name.c_str();
					cout << name <<endl;
					FileWriter::DotFileWriter(digraphs, numTask, p);

					delete[] digraphs;
				}
				catch(bad_alloc& ba)
				{
					cerr << "Bad allocation!" << endl;
					exit(EXIT_FAILURE);
					delete[] digraphs;
				}
			}
		}
	}
}

void DigraphExperiments::schedAnalysis(DigraphMethodChoice mc, string file, string directory,int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun) {
	// write a file
	string result = file+digraphFuncNames[mc];
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
	Timer timer;

	for (int num=minNum; num<=maxNum; num+=stepNum) {
		for (int util = minUtil; util<=maxUtil; util+=stepUtil) {
			int nException = 0;
			int nUnSuccess = 0;

			vector<double> recTimes;
			double sumTUtil2 = 0;

			for (int run=startRun; run<=endRun; run++) {
				string input = directory+"\\digraphNum"+Utility::int_to_string(num)+"Util"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run)+".dot";
				const char *p = input.c_str();
				cout << input <<endl;
				fout << "Read file " << input << endl;
				Digraph** digraphs;
				int taskNum;
				FileReader::DotFileReader(digraphs,taskNum, 1, p);

				try {
					timer.start();
					bool schedulable = doDigraphSchedAnalysis(mc,digraphs,num);
					timer.end();

					cout<<file<<"=>"<< DigraphExperiments::digraphFuncNames[mc] << "\t" << schedulable << endl;
					fout<<file<<"=>"<< DigraphExperiments::digraphFuncNames[mc] << "\t" << schedulable << endl;

					double tUtil = 0;
					for (int i=0; i<num; i++) {
						Digraph* digraph = digraphs[i];
						tUtil+=digraph->linear_factor;
					}


					cout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil <<endl;
					fout<<"Expected utilization = "<<1.0*util/100<<", Actual utilization = "<<tUtil <<endl;

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;
				} 
				catch(bad_alloc& ba)
				{
					nException++;
					cout << "================bad allocation===================" << endl;

					exit(EXIT_FAILURE);

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;
				}
			}


			fout<<"nException="<<nException<<endl;

			SchedulabilityAnalysis::nStateflows = num;
			SchedulabilityAnalysis::totalUtilization = util;
			// reset
			SchedulabilityAnalysis::reset();

			fout<<"Chocice="<<choice<<endl;

			SchedulabilityAnalysis::output_vectors(fout);

			fout << "Output the sum of the nodes, edges, strongly connected digraphs." << endl;
			int nNode = 0, nEdge = 0, nSSD=0;
			for (int i=0; i<SchedulabilityAnalysis::vec_nDigraphEdge.size(); i++) {
				nNode += SchedulabilityAnalysis::vec_nDigraphNode.at(i);
				nEdge += SchedulabilityAnalysis::vec_nDigraphEdge.at(i);
				nSSD += SchedulabilityAnalysis::vec_nIrreducibleDigraphs.at(i);
			}

			fout << nNode << "\t" << nEdge << "\t" << nSSD << endl;

			fout << "=============================================" << endl;
		}
	}
	fout<<"Chocice="<<choice<<endl;

	SchedulabilityAnalysis::output_vectors(fout);

	fout << "Output the sum of the nodes, edges, strongly connected digraphs." << endl;
	int nNode = 0, nEdge = 0, nSSD=0;
	for (int i=0; i<SchedulabilityAnalysis::vec_nDigraphEdge.size(); i++) {
		nNode += SchedulabilityAnalysis::vec_nDigraphNode.at(i);
		nEdge += SchedulabilityAnalysis::vec_nDigraphEdge.at(i);
		nSSD += SchedulabilityAnalysis::vec_nIrreducibleDigraphs.at(i);
	}

	fout << nNode << "\t" << nEdge << "\t" << nSSD << endl;

	fout.close();
}	

void DigraphExperiments::generateDigraphsForTestingImprovementOnTF(string directory, int minUtil, int maxUtil, int stepUtil, int totalRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int numTask;
		for (int util = minUtil; util <= maxUtil; util += stepUtil) {
			int index = 0;
			while (index < totalRun) {
				Digraph** digraphs = RandomGenerator::generate_digraph_system4(numTask,0.01*util);

				double realUtil = 0;
				for (int i=0; i<numTask; i++) {
					Digraph* digraph = digraphs[i];
					realUtil += digraph->linear_factor;
				}

				cout << "ExpectedUtil=" << 0.01*util << "\tRealUtil=" << realUtil << endl;

				if (abs(0.01*util-realUtil) <= 0.01) {
					string name = directory+"\\digraphUtil"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(index)+".dot";
					const char *p = name.c_str();
					cout << name <<endl;
					cout << numTask << endl;
					FileWriter::DotFileWriter(digraphs, numTask, p);
					index++;
				}

				//digraphs[0]->write_graphviz(cout);

				for (int i=0; i<numTask; i++) {
					//cout << "Relese digraph " << i <<endl;
					delete digraphs[i];
				}
				delete[] digraphs;
			}
		}
}

void DigraphExperiments::testImprovementOnTF(string directory, string file, int minUtil, int maxUtil, int stepUtil, int totalRun) {
	const char* p = file.c_str();

	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	SchedulabilityAnalysis::set_zero();

	vector<double> vec_sum_tf1;
	vector<double> vec_sum_tf2;

	for (int util = minUtil; util <= maxUtil; util += stepUtil) {
		double sum_tf1 = 0;
		double sum_tf2 = 0;

		for (int run=0; run<totalRun; run++) {
			//string input = directory+"\\digraphUtil"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run);
			string input = directory+"\\digraphUtil"+Utility::int_to_string(util)+"Run"+Utility::int_to_string(run)+".dot";
			const char *p = input.c_str();
			cout << input <<endl;
			fout << "Read file " << input << endl;
			Digraph** digraphs;
			int num;
			FileReader::DotFileReader(digraphs,num, 1, p);

			SchedulabilityAnalysis::prepare_all_digraphs2(digraphs,num);

			for (int i=0; i<num; i++) {
				Digraph* digraph = digraphs[i];
				if (i==num-1) {
					cout << "Improvement=" << (digraph->tf1-digraph->tf2)/digraph->tf1 << endl;
					fout << "Improvement=" << (digraph->tf1-digraph->tf2)/digraph->tf1 << endl;
					sum_tf1 += digraph->tf1;
					sum_tf2 += digraph->tf2;
				}
			}

			for (int i=0; i<num; i++)
				delete digraphs[i];
			delete[] digraphs;
		}

		vec_sum_tf1.push_back(sum_tf1);
		vec_sum_tf2.push_back(sum_tf2);
	}

	fout << "sumTF1=[";
	for (vector<double>::iterator iter = vec_sum_tf1.begin(); iter != vec_sum_tf1.end(); iter ++) {
		fout << *iter << ",";
	}
	fout << "]" << endl;

	fout << "sumTF2=[";
	for (vector<double>::iterator iter = vec_sum_tf2.begin(); iter != vec_sum_tf2.end(); iter ++) {
		fout << *iter << ",";
	}
	fout << "]" << endl;

	fout.close();
}

void DigraphExperiments::generateOneDigraphForTestingImprovemntOnIBFCalculation(string directory,int periodType, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int exception = 0;
	for (int run=0; run<maxRun; run++) {
		Digraph* digraph;
		try {
			int times = 0;
			while(1) {
				digraph = RandomGenerator::generate_DIGRAPH_VII(0,periodType);
				digraph->generate_strongly_connected_components();
				digraph->check_strongly_connected();

				if (digraph->strongly_connected) break;
				times++;
				cout << "times=" << times << endl;

				delete digraph;
			}

			string name = directory + "\\digraphType" +Utility::int_to_string(periodType)+"Run"+ Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			cout << name <<endl;
			FileWriter::DotFileWriter(digraph, p);

			delete digraph;
		}
		catch (bad_alloc) {
			exception++;
			cout << "bad_alloc exception. " << exception << endl;
			delete digraph;
		}
	}
}

void DigraphExperiments::generateOneDigraphForTestingImprovemntOnIBFCalculation2(string directory, int utilInteger, int numNode, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	double util = 1.0*utilInteger/100;

	int exception = 0;
	for (int run=0; run<maxRun; run++) {
		Digraph* digraph;
		try {
			int times = 0;
			while(1) {
				digraph = RandomGenerator::generate_one_digraph(0,10,numNode);

				int base = RandomGenerator::calculate_base();

				// set random separation times for the edges
				for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
					Edge* edge = *iter;
					int separationTime = RandomGenerator::calculate_separation_time(base);
					edge->set_separation_time(separationTime*10);
				}

				// set deadlines and wcets for the nodes
				for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
					Node* node = *iter;

					int minSeparationTime = INT_MAX;
					for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) {
						minSeparationTime = min(minSeparationTime, (*e)->separationTime);
					}

					node->set_deadline(minSeparationTime);

					int wcet = rand()%(minSeparationTime)+1;
					//int wcet = util*minSeparationTime;
					//if (wcet == 0) wcet = 1;
					node->set_wcet(wcet);
				}

				/// scale wcets to the expected utilization
				// step 1: calculate the linear factor
				digraph->generate_strongly_connected_components();
				digraph->check_strongly_connected();

				if (!digraph->strongly_connected) {
					delete digraph;
					times ++;
					cout << "times=" << times << endl;
					continue;
				}
#if 0 // randomly set the wcets
				digraph->calculate_period_gcd();
				digraph->calculate_linear_factor();
				//cout<<"linear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
				// step 2: rescale wcets
				double factor = util/digraph->linear_factor;
				digraph->scale_wcet(factor);
#endif
				break;
			}

			string name = directory + "\\digraphUtil" +Utility::int_to_string(utilInteger)+"NodeNum"+Utility::int_to_string(numNode)+"Run"+ Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			cout << name <<endl;
			FileWriter::DotFileWriter(digraph, p);

			delete digraph;
		}
		catch (bad_alloc) {
			exception++;
			cout << "bad_alloc exception. " << exception << endl;
			delete digraph;
		}
	}
}

void DigraphExperiments::generateOneDigraphForTestingImprovemntOnIBFCalculation2(string directory, int utilInteger, int minNode, int maxNode, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	double util = 1.0*utilInteger/100;

	int exception = 0;
	for (int run=0; run<maxRun; run++) {
		Digraph* digraph;
		int numNode = minNode + rand()% (maxNode-numNode+1);
		cout << "numNode = " << numNode << endl;
		try {
			int times = 0;
			while(1) {
				digraph = RandomGenerator::generate_one_digraph(0,10,numNode);

				int base = RandomGenerator::calculate_base();

				// set random separation times for the edges
				for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
					Edge* edge = *iter;
					int separationTime = RandomGenerator::calculate_separation_time(base);
					edge->set_separation_time(separationTime*10);
				}

				// set deadlines and wcets for the nodes
				for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
					Node* node = *iter;

					int minSeparationTime = INT_MAX;
					for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) {
						minSeparationTime = min(minSeparationTime, (*e)->separationTime);
					}

					node->set_deadline(minSeparationTime);

					int wcet = rand()%minSeparationTime+1;
					//int wcet = util*minSeparationTime;
					//if (wcet == 0) wcet = 1;
					node->set_wcet(wcet);
				}

				/// scale wcets to the expected utilization
				// step 1: calculate the linear factor
				digraph->generate_strongly_connected_components();
				digraph->check_strongly_connected();

				if (!digraph->strongly_connected) {
					delete digraph;
					times ++;
					cout << "times=" << times << endl;
					continue;
				}
#if 0 // randomly set the wcets
				digraph->calculate_period_gcd();
				digraph->calculate_linear_factor();
				//cout<<"linear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
				// step 2: rescale wcets
				double factor = util/digraph->linear_factor;
				digraph->scale_wcet(factor);
#endif
				break;
			}

			string name = directory + "\\digraphUtil" +Utility::int_to_string(utilInteger)+"NodeNum"+Utility::int_to_string(1025)+"Run"+ Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			cout << name <<endl;
			FileWriter::DotFileWriter(digraph, p);

			delete digraph;
		}
		catch (bad_alloc) {
			exception++;
			cout << "bad_alloc exception. " << exception << endl;
			delete digraph;
		}
	}
}

void DigraphExperiments::testImprovemntOnIBFCalculationInOneDigrah(string directory, string file, int periodType, int minRun, int maxRun, int stepRun) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int exception = 0;
	Timer timer;
	Digraph* digraph;

	try {
		map<int,double> map_sumTCalIBFWithPeriodicity, map_sumTCalIBFWithoutPeriodicity;
		map<int,double> map_sumTCalRBFWithPeriodicity, map_sumTCalRBFWithoutPeriodicity;
		map<int,int> map_sumSuccess;
		int sumGDRTSize = 0;
		int sumUDRTSize = 0;

		for (int run=minRun; run<=maxRun; run+=stepRun) {

			string name = directory + "\\digraphType" +Utility::int_to_string(periodType)+"Run"+ Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			FileReader::DotFileReader(digraph,1,p);
			cout << name << endl;

			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			if (!digraph->strongly_connected) {
				cerr << "The generated digraph has to be strongly connected." << endl;
				exit(EXIT_FAILURE);
			}

			digraph->calculate_period_gcd();
			digraph->calculate_all_gcd();
			digraph->calculate_linear_factor();
			//digraph->tf = 50000; // tf is large enough to calculate linear defect

			cout << "lfac=" << digraph->linear_factor << endl;

			map<int,double> map_tCalIBFWithPeriodicity, map_tCalIBFWithoutPeriodicity;
			map<int,double> map_tCalRBFWithPeriodicity, map_tCalRBFWithoutPeriodicity;
			map<int,int> map_success;

			double tCalIBFLinearization;
			double tCalIBFNonLinearization;
			double tCalRBFLinearization;
			double tCalRBFNonLinearization;
			bool hasCalculatedIBFLinearization = false;
			bool hasCalculatedIBFNonLinearization = false;
			bool hasCalculatedRBFLinearization = false;
			bool hasCalculatedRBFNonLinearization = false;
			map<int,double> calIBFTimes;
			map<int,double> calRBFTimes;

			int startTime = 5000;
			int endTime = 100000;
			int stepTime = 5000;


			for (int tf=startTime; tf<=endTime; tf += stepTime) {
			//for (int tf=5000; tf<=5000; tf += 5000) {
				// calculate IBF with periodicity
				cout << "Calculating ibf with periodicity..." << endl;
				digraph->tf = 10000;
				map<int,double> ibf_map;
				map<int,int> rbf_map;

				bool success = digraph->calculate_ibf_with_periodicity(false,hasCalculatedIBFLinearization,tf, ibf_map, tCalIBFLinearization, tCalIBFNonLinearization);
				digraph->calculate_rbf_with_periodicity(false,hasCalculatedRBFLinearization,tf, rbf_map, tCalRBFLinearization, tCalRBFNonLinearization);
				//digraph->prepare_ibf_calculation(false);
				double tCalIBFWithPeriodicity = tCalIBFLinearization + tCalIBFNonLinearization;
				double tCalRBFWithPeriodicity = tCalRBFLinearization + tCalRBFNonLinearization;
				map_tCalIBFWithPeriodicity[tf] = tCalIBFWithPeriodicity;
				map_tCalRBFWithPeriodicity[tf] = tCalRBFWithPeriodicity;

				if (map_sumTCalIBFWithPeriodicity.find(tf) == map_sumTCalIBFWithPeriodicity.end()) {
					map_sumTCalIBFWithPeriodicity[tf] = tCalIBFWithPeriodicity;
					map_sumTCalRBFWithPeriodicity[tf] = tCalRBFWithPeriodicity;
				}
				else {
					map_sumTCalIBFWithPeriodicity[tf] += tCalIBFWithPeriodicity;
					map_sumTCalRBFWithPeriodicity[tf] += tCalRBFWithPeriodicity;
				}

				cout << "tf=" << tf << "\t" << "tCalIBFWithPeriodicity=" << tCalIBFWithPeriodicity << "tCalRBFWithPeriodicity=" << tCalRBFWithPeriodicity << endl;

				if (success) map_success[tf] = 1;
				else map_success[tf] = 0;

				if (map_sumSuccess.find(tf) == map_sumSuccess.end()) map_sumSuccess[tf] = map_success[tf];
				else map_sumSuccess[tf] += map_success[tf];


				cout << "Calculating ibf without periodicity..." << endl;
				map<int, int> ibf_map2;
				map<int, int> rbf_map2;
				timer.start();
				digraph->calculate_ibf_without_periodicity2(tf, ibf_map2);
				timer.end();
				//digraph->calculate_ibf_without_periodicity(hasCalculatedNonLinearization, startTime, endTime, stepTime, ibf_map2, calTimes);
				double tCalIBFWithoutPeriodicity = timer.getTime(); // calTimes[tf];
				map_tCalIBFWithoutPeriodicity[tf] = tCalIBFWithoutPeriodicity;

				timer.start();
				digraph->calculate_rbf_without_periodicity2(tf, rbf_map2);
				timer.end();
				//digraph->calculate_ibf_without_periodicity(hasCalculatedNonLinearization, startTime, endTime, stepTime, ibf_map2, calTimes);
				double tCalRBFWithoutPeriodicity = timer.getTime(); // calTimes[tf];
				map_tCalRBFWithoutPeriodicity[tf] = tCalRBFWithoutPeriodicity;

				if (map_sumTCalIBFWithoutPeriodicity.find(tf) == map_sumTCalIBFWithoutPeriodicity.end()) {
					map_sumTCalIBFWithoutPeriodicity[tf] = tCalIBFWithoutPeriodicity;
					map_sumTCalRBFWithoutPeriodicity[tf] = tCalRBFWithoutPeriodicity;
				}
				else {
					map_sumTCalIBFWithoutPeriodicity[tf] += tCalIBFWithoutPeriodicity;
					map_sumTCalRBFWithoutPeriodicity[tf] += tCalRBFWithoutPeriodicity;
				}

				cout << "tf=" << tf << "\t" << "tCalIBFWithoutPeriodicity=" << tCalIBFWithoutPeriodicity << "tCalRBFWithoutPeriodicity=" << tCalRBFWithoutPeriodicity  << endl;
			}

			sumGDRTSize += digraph->gran_digraph->n_size;
			sumUDRTSize += digraph->unit_digraph->n_size;

			cout << "sumGDRTSize=" << sumGDRTSize << "\tsumUDRTSize=" << sumUDRTSize <<endl;
			fout << "sumGDRTSize=" << sumGDRTSize << "\tsumUDRTSize=" << sumUDRTSize <<endl;
			fout << "Output the IBF calculating time with periodicity..." << endl;
			for (map<int,double>::iterator iter = map_tCalIBFWithPeriodicity.begin(); iter != map_tCalIBFWithPeriodicity.end(); iter++) {
				fout << "tf=" << iter->first << "\t" << "tCalIBFWithPeriodicity=" << iter->second << endl;
			}

			fout << "Output the IBF calculating time without periodicity..." << endl;
			for (map<int,double>::iterator iter = map_tCalIBFWithoutPeriodicity.begin(); iter != map_tCalIBFWithoutPeriodicity.end(); iter++) {
				fout << "tf=" << iter->first << "\t" << "tCalIBFWithoutPeriodicity=" << iter->second << endl;
			}

			fout << "Output the RBF calculating time with periodicity..." << endl;
			for (map<int,double>::iterator iter = map_tCalRBFWithPeriodicity.begin(); iter != map_tCalRBFWithPeriodicity.end(); iter++) {
				fout << "tf=" << iter->first << "\t" << "tCalRBFWithPeriodicity=" << iter->second << endl;
			}

			fout << "Output the RBF calculating time without periodicity..." << endl;
			for (map<int,double>::iterator iter = map_tCalRBFWithoutPeriodicity.begin(); iter != map_tCalRBFWithoutPeriodicity.end(); iter++) {
				fout << "tf=" << iter->first << "\t" << "tCalRBFWithoutPeriodicity=" << iter->second << endl;
			}

			cout << "GranDigraph=>\tlfac=" << digraph->gran_digraph->lfac 
				<< "\tlper=" << digraph->gran_digraph->lper 
				<< "\tldef=" << digraph->gran_digraph->ldef
				<< "\tubldef=" << digraph->gran_digraph->ubldef << endl;

			fout << "====================" << run+1 << "==============" << endl;
			cout << "All: output the IBF calculating time with periodicity..." << endl;
			fout << "sumTCalIBFWithPeriodicity=[";
			for (map<int,double>::iterator iter = map_sumTCalIBFWithPeriodicity.begin(); iter != map_sumTCalIBFWithPeriodicity.end(); iter++) {
				cout << "tf=" << iter->first << "\t" << "sumTCalIBFWithPeriodicity=" << iter->second << endl;
				fout << iter->second << ",";
			}
			fout << "]" << endl;

			cout << "All: output the IBF calculating time without periodicity..." << endl;
			fout << "sumTCalIBFWithoutPeriodicity=[";
			for (map<int,double>::iterator iter = map_sumTCalIBFWithoutPeriodicity.begin(); iter != map_sumTCalIBFWithoutPeriodicity.end(); iter++) {
				cout << "tf=" << iter->first << "\t" << "sumTCalIBFWithoutPeriodicity=" << iter->second << endl;
				fout << iter->second << ",";
			}
			fout << "]" << endl;

			cout << "All: output the RBF calculating time with periodicity..." << endl;
			fout << "sumTCalRBFWithPeriodicity=[";
			for (map<int,double>::iterator iter = map_sumTCalRBFWithPeriodicity.begin(); iter != map_sumTCalRBFWithPeriodicity.end(); iter++) {
				cout << "tf=" << iter->first << "\t" << "sumTCalRBFWithPeriodicity=" << iter->second << endl;
				fout << iter->second;
				if (iter != (--map_sumTCalRBFWithPeriodicity.end()))
					fout << ",";
			}
			fout << "]" << endl;

			cout << "All: output the RBF calculating time without periodicity..." << endl;
			fout << "sumTCalRBFWithoutPeriodicity=[";
			for (map<int,double>::iterator iter = map_sumTCalRBFWithoutPeriodicity.begin(); iter != map_sumTCalRBFWithoutPeriodicity.end(); iter++) {
				cout << "tf=" << iter->first << "\t" << "sumTCalRBFWithoutPeriodicity=" << iter->second << endl;
				fout << iter->second;
				if (iter != (--map_sumTCalRBFWithoutPeriodicity.end()))
					fout << ",";
			}
			fout << "]" << endl;

			delete digraph;
		}

		fout << "====================ALL====================" << endl;
		fout << "sumGDRTSize=" << sumGDRTSize << "\tsumUDRTSize=" << sumUDRTSize <<endl;
		cout << "All: output the IBF calculating time with periodicity..." << endl;
		fout << "sumTCalIBFWithPeriodicity=[";
		for (map<int,double>::iterator iter = map_sumTCalIBFWithPeriodicity.begin(); iter != map_sumTCalIBFWithPeriodicity.end(); iter++) {
			cout << "tf=" << iter->first << "\t" << "sumTCalIBFWithPeriodicity=" << iter->second << endl;
			fout << iter->second;
			if (iter != (--map_sumTCalIBFWithPeriodicity.end()))
				fout << ",";
		}
		fout << "]" << endl;

		cout << "All: output the IBF calculating time without periodicity..." << endl;
		fout << "sumTCalIBFWithoutPeriodicity=[";
		for (map<int,double>::iterator iter = map_sumTCalIBFWithoutPeriodicity.begin(); iter != map_sumTCalIBFWithoutPeriodicity.end(); iter++) {
			cout << "tf=" << iter->first << "\t" << "sumTCalIBFWithoutPeriodicity=" << iter->second << endl;
			fout << iter->second;
			if (iter != (--map_sumTCalIBFWithoutPeriodicity.end()))
				fout << ",";
		}
		fout << "]" << endl;

		cout << "All: output the RBF calculating time with periodicity..." << endl;
		fout << "sumTCalRBFWithPeriodicity=[";
		for (map<int,double>::iterator iter = map_sumTCalRBFWithPeriodicity.begin(); iter != map_sumTCalRBFWithPeriodicity.end(); iter++) {
			cout << "tf=" << iter->first << "\t" << "sumTCalRBFWithPeriodicity=" << iter->second << endl;
			fout << iter->second;
			if (iter != (--map_sumTCalRBFWithPeriodicity.end()))
				fout << ",";
		}
		fout << "]" << endl;

		cout << "All: output the RBF calculating time without periodicity..." << endl;
		fout << "sumTCalRBFWithoutPeriodicity=[";
		for (map<int,double>::iterator iter = map_sumTCalRBFWithoutPeriodicity.begin(); iter != map_sumTCalRBFWithoutPeriodicity.end(); iter++) {
			cout << "tf=" << iter->first << "\t" << "sumTCalRBFWithoutPeriodicity=" << iter->second << endl;
			fout << iter->second;
			if (iter != (--map_sumTCalRBFWithoutPeriodicity.end()))
				fout << ",";
		}
		fout << "]" << endl;

		cout << "All: success..." << endl;
		fout << "sumSuccess=[";
		for (map<int,int>::iterator iter = map_sumSuccess.begin(); iter != map_sumSuccess.end(); iter++) {
			cout << "tf=" << iter->first << "\t" << "sumSuccess=" << iter->second << endl;
		}
		fout << "]" << endl;
	}
	catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
		delete digraph;
	}

	fout.close();
}

void DigraphExperiments::testImprovemntOnIBFCalculationInOneDigrah2(string directory, string file, int util, int numNode, int minRun, int maxRun, int stepRun) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int exception = 0;
	Timer timer;
	Digraph* digraph;

	vector<double> ratios;
	vector<double> ratios2;
	vector<double> tPs;
	vector<double> tNPs;

	try {
		for (int nFactor=10; nFactor<=100; nFactor += 10) {
			double sum_ratio = 0;
			double sum_t0 = 0;
			double sum_t1 = 0;
			for (int run=minRun; run<=maxRun; run+=stepRun) {
				string name = directory + "\\digraphUtil" +Utility::int_to_string(util)+"NodeNum"+Utility::int_to_string(numNode)+"Run"+ Utility::int_to_string(run) + ".dot";
				const char *p = name.c_str();
				FileReader::DotFileReader(digraph,1,p);
				cout << name << endl;
				fout << name << endl;

				double tP = 0;
				double tNP = 0;

				timer.start();
				digraph->generate_strongly_connected_components();
				digraph->check_strongly_connected();
				timer.end();
				tP += timer.getTime();

				if (!digraph->strongly_connected) {
					cerr << "The generated digraph has to be strongly connected." << endl;
					exit(EXIT_FAILURE);
				}

				int maxSeperationTime = INT_MIN;
				for (auto edge : digraph->edge_vec) {
					maxSeperationTime = max(maxSeperationTime,edge->separationTime);
				}
				int tf = nFactor*maxSeperationTime;
				cout << "tf = " << tf << endl;
				fout << "tf = " << tf << endl;

				timer.start();
				digraph->calculate_period_gcd();
				timer.end();
				tP += timer.getTime();
				tNP += timer.getTime();

				timer.start();
				digraph->calculate_all_gcd();
				digraph->calculate_linear_factor();
				timer.end();
				tP += timer.getTime();
				//digraph->tf = 50000; // tf is large enough to calculate linear defect

				cout << "lfac=" << digraph->linear_factor << endl;
				fout << "lfac=" << digraph->linear_factor << endl;

				map<int,int> ibf_map0;
				map<int,int> ibf_map1;
			
				// calculate IBF with periodicity
				cout << "Calculating ibf with periodicity..." << endl;
				timer.start();
				digraph->unit_digraph = new UnitDigraph(digraph);
				digraph->unit_digraph->prepare3(false);
				digraph->calculate_ibf_with_periodicity_DP(tf,ibf_map0);
				timer.end();
				tP += timer.getTime();

				cout << "Calculating ibf without periodicity..." << endl;
				timer.start();
				digraph->calculate_ibf_without_periodicity_DP(tf,ibf_map1);	
				timer.end();
				tNP += timer.getTime();

				cout << "nFactor = " << nFactor << ", run = " << run << ", ratio = " << tNP/tP << ", tP = " << tP << ", tNP = " << tNP  << endl; 
				fout << "nFactor = " << nFactor << ", run = " << run << ", ratio = " << tNP/tP << ", tP = " << tP << ", tNP = " << tNP  << endl;

				sum_ratio += tNP/tP;

				sum_t0 += tP;
				sum_t1 += tNP;

				delete digraph;
			}
			ratios.push_back(sum_ratio/(maxRun-minRun+1));
			ratios2.push_back(sum_t1/sum_t0);
			tPs.push_back(sum_t0);
			tNPs.push_back(sum_t1);
		}
	}
	catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
		delete digraph;
	}

	Utility::output_one_vector(cout,"ratios",ratios);
	Utility::output_one_vector(fout,"ratios",ratios);

	Utility::output_one_vector(cout,"ratios2",ratios2);
	Utility::output_one_vector(fout,"ratios2",ratios2);

	Utility::output_one_vector(cout,"tPs",tPs);
	Utility::output_one_vector(fout,"tPs",tPs);

	Utility::output_one_vector(cout,"tNPs",tNPs);
	Utility::output_one_vector(fout,"tNPs",tNPs);

	fout.close();
}

void DigraphExperiments::testLinearDefectForIBFRBF(string directory, string file, int periodType, int minRun, int maxRun, int stepRun) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}
	
	int exception = 0;
	Timer timer;
	Digraph* digraph;

	try {
		for (int run=minRun; run<=maxRun; run+=stepRun) {

			string name = directory + "\\digraphType" +Utility::int_to_string(periodType)+"Run"+ Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			FileReader::DotFileReader(digraph,1,p);
			cout << name << endl;

			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			if (!digraph->strongly_connected) {
				cerr << "The generated digraph has to be strongly connected." << endl;
				exit(EXIT_FAILURE);
			}

			digraph->calculate_period_gcd();
			digraph->calculate_all_gcd();
			digraph->calculate_linear_factor();
			digraph->tf = 50000; // tf is large enough to calculate linear defect

			cout << "lfac=" << digraph->linear_factor << endl;

			digraph->prepare_ibf_calculation(false);
			digraph->prepare_rbf_calculation(false);

			cout << name << "=>ldef of ibf = " << digraph->gran_digraph->ldef << "\t ldef of rbf=" << digraph->unit_digraph->ldef << endl;
			fout << name << "=>ldef of ibf = " << digraph->gran_digraph->ldef << "\t ldef of rbf=" << digraph->unit_digraph->ldef << endl;

			delete digraph;
		}
	}
	catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
		delete digraph;
	}

	fout.close();
}

void DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriod(string directory, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int exception = 0;

	Digraph* hDigraph;

	try {
		for (int run=0; run<maxRun; run++) {
			// first we generate one digraph with around 0.04 utilization

			int times = 0;
			while(1) {
				hDigraph = RandomGenerator::generate_DIGRAPH_IV(0,0.04);
				hDigraph->generate_strongly_connected_components();
				hDigraph->check_strongly_connected();

				if (!hDigraph->strongly_connected) {
					times++;
					cout << "Not strongly connected times=" << times << endl;
					delete hDigraph;

					continue;
				}

				hDigraph->calculate_linear_factor();
				cout << "\tExpectedUtil=" << 0.04 << "\tRealUtil=" << hDigraph->linear_factor << endl;
				
				if (abs(0.04-hDigraph->linear_factor) < 0.003) {
					break; // we find a digraph satisfied the need
				}

				delete hDigraph;
			}

			string name = directory + "\\digraph" + Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			cout << name <<endl;
			FileWriter::DotFileWriter(hDigraph, p);

			delete hDigraph;
		}

	} catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
		delete hDigraph;
	}

}

void DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriod2(string directory, int fixedNodeNum, int minRun, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int exception = 0;

	Digraph* hDigraph;

	try {
		for (int run=minRun; run<=maxRun; run++) {
			// first we generate one digraph with around 0.04 utilization
			int times = 0;
			while (1) {
				hDigraph = RandomGenerator::generate_one_digraph(0,RandomGenerator::DIGRAPH_SCALE,fixedNodeNum,0.04);
				//hDigraph->generate_strongly_connected_components();
				//hDigraph->check_strongly_connected();

				if (hDigraph->strongly_connected) break;

				cout << "Run" << ++times << endl;
				delete hDigraph;
			}

			string name = directory + "\\digraph" + Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			cout << name <<endl;
			FileWriter::DotFileWriter(hDigraph, p);

			delete hDigraph;
		}

	} catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
		delete hDigraph;
	}
}

void DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriodUtil(string directory, int fixedNodeNum, int minUtil, int maxUtil,int stepUtil, int minRun, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int exception = 0;

	Digraph* hDigraph;

	try {
		for (int util = minUtil; util <=maxUtil; util += stepUtil) {
			for (int run=minRun; run<=maxRun; run++) {
				// first we generate one digraph with around 0.04 utilization
				int times = 0;
				while (1) {
					hDigraph = RandomGenerator::generate_one_digraph(0,RandomGenerator::DIGRAPH_SCALE,fixedNodeNum,0.01*util/5);
					//hDigraph->generate_strongly_connected_components();
					//hDigraph->check_strongly_connected();

					if (hDigraph->strongly_connected) break;

					cout << "Run" << ++times << endl;
					delete hDigraph;
				}

				string name = directory + "\\digraphUtil" + Utility::int_to_string(util) + "Run" + Utility::int_to_string(run) + ".dot";
				const char *p = name.c_str();
				cout << name <<endl;
				FileWriter::DotFileWriter(hDigraph, p);

				delete hDigraph;
			}
		}
	} catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
		delete hDigraph;
	}
}

void DigraphExperiments::testEffectOfDigraphPeriod(string directory, string file, int startRun, int endRun, int choice) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int exception = 0;

	vector<vector<double>> AllIBFResults;
	vector<vector<double>> AllRBFResults;
	vector<vector<double>> AllLUBRBFResults;
	vector<vector<double>> AllLUBIBFResults;

	try {
		int numNonSucc = 0;
		for (int run=startRun; run<=endRun; run++) {
			vector<vector<double>> IBFResults;
			vector<vector<double>> RBFResults;
			vector<vector<double>> LUBRBFResults;
			vector<vector<double>> LUBIBFResults;

			string name = directory + "\\digraph" + Utility::int_to_string(run) + ".dot";
			const char *p = name.c_str();
			cout << name << endl;

			bool success = true;

			// [0.2,0.4,0.6,0.8]
			for (int util = 2; util <=8; util+=2) {
			//for (int util = 10; util <=10; util+=2) {
				vector<double> IBFResult;
				vector<double> RBFResult;
				vector<double> LUBRBFResult;
				vector<double> LUBIBFResult;

				double exactUtil;
				set<double> ratios;

				if (choice == 1) 
					success = testEffectOfDigraphPeriod(IBFResult,RBFResult,LUBRBFResult,LUBIBFResult,ratios, p,util,exactUtil);
				else if (choice == 2)
					testEffectOfDigraphPeriod2(IBFResult,LUBRBFResult,LUBIBFResult, ratios,p,util,exactUtil);
				else if (choice == 3)
					success = testEffectOfDigraphPeriodWithHarmonicSet(IBFResult,RBFResult,LUBRBFResult,LUBIBFResult,ratios, p,util,exactUtil);
				else {
					cerr << "Error choice for testEffectOfDigraphPeriod, choice=" << choice << endl;
					exit(EXIT_FAILURE);
				}

				if (!success) {
					numNonSucc++;
					break;
				}

				cout << "run=" << run << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << exactUtil << endl;

				fout << "the number of ratios = " << ratios.size() << endl;
				fout << "ratios=[";
				for (set<double>::iterator iter = ratios.begin(); iter != ratios.end(); iter++) {
					fout << *iter;  
					if (iter != (--ratios.end()))
						fout << ",";
				}
				fout << "];" << endl;

				IBFResults.push_back(IBFResult);
				RBFResults.push_back(RBFResult);
				LUBRBFResults.push_back(LUBRBFResult);
				LUBIBFResults.push_back(LUBIBFResult);
			}

			if (!success) continue;

			/*
			// output
			outputResults(fout, IBFResults, "ratioExact");
			outputResults(fout, LUBRBFResults, "ratioLRBF");
			outputResults(fout, LUBIBFResults, "ratioLIBF");
			*/

			if (AllIBFResults.empty()) {
				AllIBFResults = IBFResults;
				AllRBFResults = RBFResults;
				AllLUBIBFResults = LUBIBFResults;
				AllLUBRBFResults = LUBRBFResults;
			}
			else {
				for ( int i=0; i<AllIBFResults.size(); i++) {
					for (int j=0; j<AllIBFResults.at(i).size(); j++) {
						AllIBFResults.at(i).at(j) = AllIBFResults.at(i).at(j) + IBFResults.at(i).at(j);
						AllRBFResults.at(i).at(j) = AllRBFResults.at(i).at(j) + RBFResults.at(i).at(j);
						AllLUBIBFResults.at(i).at(j) = AllLUBIBFResults.at(i).at(j) + LUBIBFResults.at(i).at(j);
						AllLUBRBFResults.at(i).at(j) = AllLUBRBFResults.at(i).at(j) + LUBRBFResults.at(i).at(j);
					}
				}
			}

			fout << "======================Run=" << run+1 << "\tNumNonSucc=" << numNonSucc << "=====================" << endl;
			outputResults(fout, AllIBFResults, "ratioIBF");
			fout << endl;
			outputResults(fout, AllRBFResults, "ratioRBF");
			fout << endl;
			outputResults(fout, AllLUBRBFResults, "ratioLRBF");
			fout << endl;
			outputResults(fout, AllLUBIBFResults, "ratioLIBF");
			fout << endl;
		}

		fout << "======================ALL=" << endRun-startRun+1 << "\tNumNonSucc=" << numNonSucc << "=====================" << endl;
		outputResults(fout, AllIBFResults, "ratioIBF");
		fout << endl;
		outputResults(fout, AllRBFResults, "ratioRBF");
		fout << endl;
		outputResults(fout, AllLUBRBFResults, "ratioLRBF");
		fout << endl;
		outputResults(fout, AllLUBIBFResults, "ratioLIBF");
		fout << endl;

	} catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
	}

}

void DigraphExperiments::testEffectOfDigraphPeriodUtil(string directory, string file, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int exception = 0;

	vector<vector<double>> AllIBFResults;
	vector<vector<double>> AllRBFResults;
	vector<vector<double>> AllLUBRBFResults;
	vector<vector<double>> AllLUBIBFResults;
	vector<int> AllNumNonSucc;

	set<double> ratios;
	ratios.insert(0.2);
	ratios.insert(0.5);
	ratios.insert(0.8);

	try {
		for (int util = minUtil; util <=maxUtil; util+=stepUtil) {
			vector<double> IBFResults;
			vector<double> RBFResults;
			vector<double> LUBRBFResults;
			vector<double> LUBIBFResults;

			bool success = true;
			int numNonSucc = 0;

			for (int run=startRun; run<=endRun; run++) {
				string name = directory + "\\digraphUtil" + Utility::int_to_string(util) + "Run" + Utility::int_to_string(run) + ".dot";
				const char *p = name.c_str();
				cout << name << endl;

				vector<double> IBFResult;
				vector<double> RBFResult;
				vector<double> LUBRBFResult;
				vector<double> LUBIBFResult;

				double exactUtil;
				
				success = testEffectOfDigraphPeriodWithHarmonicSetUtil(IBFResult,RBFResult,LUBRBFResult,LUBIBFResult,ratios, p,util,exactUtil);

				if (!success) {
					numNonSucc++;
					continue;
				}

				cout << "run=" << run << "\tExpectedUtil=" << 0.01*util/5 << "\tExactUtil=" << exactUtil << endl;

				if (IBFResults.empty()) {
					IBFResults = IBFResult;
					RBFResults = RBFResult;
					LUBRBFResults = LUBRBFResult;
					LUBIBFResults = LUBIBFResult;
				} else {
					for (int i=0; i<IBFResults.size(); i++) {
						IBFResults[i] = IBFResults[i] + IBFResult[i];
						RBFResults[i] = RBFResults[i] + RBFResult[i];
						LUBIBFResults[i] = LUBIBFResults[i] + LUBIBFResult[i];
						LUBRBFResults[i] = LUBRBFResults[i] + LUBRBFResult[i];
					}
				}
			}

			AllNumNonSucc.push_back(numNonSucc);

			AllRBFResults.push_back(RBFResults);
			AllIBFResults.push_back(IBFResults);
			AllLUBIBFResults.push_back(LUBIBFResults);
			AllLUBRBFResults.push_back(LUBRBFResults);

			fout << "======================Util=" << util << "\tNumNonSucc=" << numNonSucc << "=====================" << endl;
			Utility::output_one_vector(fout, "NumNonSucc", AllNumNonSucc);
			fout << endl;
			outputResults2(fout, AllIBFResults, "ratioIBF");
			fout << endl;
			outputResults2(fout, AllRBFResults, "ratioRBF");
			fout << endl;
			outputResults2(fout, AllLUBRBFResults, "ratioLRBF");
			fout << endl;
			outputResults2(fout, AllLUBIBFResults, "ratioLIBF");
			fout << endl;
		}

		Utility::output_one_vector(fout, "NumNonSucc", AllNumNonSucc);
		fout << endl;
		outputResults2(fout, AllIBFResults, "ratioIBF");
		fout << endl;
		outputResults2(fout, AllRBFResults, "ratioRBF");
		fout << endl;
		outputResults2(fout, AllLUBRBFResults, "ratioLRBF");
		fout << endl;
		outputResults2(fout, AllLUBIBFResults, "ratioLIBF");
		fout << endl;

	} catch (bad_alloc) {
		exception++;
		cout << "bad_alloc exception. " << exception << endl;
	}

}

bool DigraphExperiments::testEffectOfDigraphPeriod(vector<double> &IBFResult, vector<double> &RBFResult, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double>& ratios,  const char *p, int util, double &exactUtil) {
	Digraph* hDigraph;
	FileReader::DotFileReader(hDigraph,1,p);

	// set the utilization
	hDigraph->generate_strongly_connected_components();
	hDigraph->check_strongly_connected();

	hDigraph->calculate_period_gcd();
	hDigraph->calculate_all_gcd();
	hDigraph->calculate_linear_factor();

	double factor = 0.02*util/hDigraph->linear_factor;
	cout << "Before scaling wcets: factor=" << factor << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << hDigraph->linear_factor << endl;

	// release sccs
	for (vector<Digraph*>::iterator iter = hDigraph->sccs.begin(); iter != hDigraph->sccs.end(); iter++) {
		(*iter)->node_vec.clear();
		(*iter)->edge_vec.clear();
		delete *iter;
		*iter = NULL;
	}
	hDigraph->sccs.clear();

	hDigraph->scale_wcet(factor);
	hDigraph->generate_strongly_connected_components();
	hDigraph->check_strongly_connected();

	hDigraph->calculate_period_gcd();
	hDigraph->calculate_all_gcd();
	hDigraph->calculate_linear_factor();

	hDigraph->calculate_linear_upper_bounds();

	cout << "After scaling wcets: factor=" << factor << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << hDigraph->linear_factor << endl;
	exactUtil = hDigraph->linear_factor;

	// calculate the periodicity parameters
	int tf_max = 500000;
	hDigraph->tf = tf_max;
	map<int,double> ibf_map;
	bool success = hDigraph->calculate_ibf_with_periodicity(false,ibf_map);

	map<int,int> rbf_map;
	hDigraph->calculate_rbf_without_periodicity(hDigraph->gran_digraph->ubldef+hDigraph->gran_digraph->lper, rbf_map);
	int lper = hDigraph->gran_digraph->lper;
	for (int t=hDigraph->gran_digraph->ubldef+hDigraph->gran_digraph->lper+1; t <= hDigraph->tf; t++) {
		rbf_map[t] = rbf_map[t-lper] + hDigraph->linear_factor*lper;  
	}
	

	if (!success) {
		cout << "We cannot calculate the periodicity parameters for ibf..." << endl;
		delete hDigraph;
		return false;
	}

	for (int i=0; i<1000; i++) {
		double ratio = 1.0-0.001*i;
		ratios.insert(ratio);
	}

	for (set<double>::iterator iter = ratios.begin(); iter != ratios.end(); iter++) {
		double ul = 0.1*util-hDigraph->linear_factor;
		double tl = 1.0*hDigraph->gran_digraph->lper/(*iter);
		double wcet = ul*tl;

		cout << "wcet = " << wcet << endl;
		cout << "lfac=" << exactUtil << "\tc_rbf = " << hDigraph->c_rbf << "\tc_ibf = " << hDigraph->c_ibf << endl;
		
		// ibf
		double rt = 0;
		for (int t=1; t<tf_max; t++) {
			rt = ibf_map[t] + wcet;
			if (rt <= t) break;
			t = rt;
			if ( t >= tf_max) {
				cout << "We have arrived at the boundary." << endl;
				delete hDigraph;
				return false;
			}
		}

		IBFResult.push_back(rt/tl);
		cout << rt << "\t";

		for (int t=1; t<tf_max; t++) {
			rt = rbf_map[t] + wcet;
			if (rt <= t) break;
			t = rt;
			if ( t >= tf_max) {
				cout << "We have arrived at the boundary." << endl;
				delete hDigraph;
				return false;
			}
		}

		RBFResult.push_back(rt/tl);
		cout << rt << "\t";

		// lubrbf
		rt = (hDigraph->c_rbf+wcet)/(1.0-hDigraph->linear_factor);
		LUBRBFResult.push_back(rt/tl);
		cout << rt << "\t";

		// lubibf
		rt = (hDigraph->c_ibf+wcet)/(1.0-hDigraph->linear_factor);
		LUBIBFResult.push_back(rt/tl);
		cout << rt << endl;
	}

	delete hDigraph;
	return true;

}

void DigraphExperiments::testEffectOfDigraphPeriod2(vector<double> &ExactResult, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double>& ratios,  const char *p, int util, double &exactUtil) {
	Digraph* hDigraph;
	FileReader::DotFileReader(hDigraph,1,p);

	// set the utilization
	hDigraph->generate_strongly_connected_components();
	hDigraph->check_strongly_connected();

	hDigraph->calculate_period_gcd();
	hDigraph->calculate_all_gcd();
	hDigraph->calculate_linear_factor();

	double factor = 0.02*util/hDigraph->linear_factor;
	cout << "Before scaling wcets: factor=" << factor << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << hDigraph->linear_factor << endl;

	// release sccs
	for (vector<Digraph*>::iterator iter = hDigraph->sccs.begin(); iter != hDigraph->sccs.end(); iter++) {
		(*iter)->node_vec.clear();
		(*iter)->edge_vec.clear();
		delete *iter;
		*iter = NULL;
	}
	hDigraph->sccs.clear();

	hDigraph->scale_wcet(factor);
	hDigraph->generate_strongly_connected_components();
	hDigraph->check_strongly_connected();

	hDigraph->calculate_period_gcd();
	hDigraph->calculate_all_gcd();
	hDigraph->calculate_linear_factor();

	hDigraph->calculate_linear_upper_bounds();

	cout << "After scaling wcets: factor=" << factor << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << hDigraph->linear_factor << endl;
	exactUtil = hDigraph->linear_factor;

	// calculate the periodicity parameters
	hDigraph->tf = 500000;
	map<int,double> ibf_map;
	bool success = hDigraph->calculate_ibf_with_periodicity(false,ibf_map);

	if (!success) {
		cerr << "We cannot calculate the periodicity parameters for ibf..." << endl;
		exit(EXIT_FAILURE);
	}


	for (int i=1; i<=10; i++) {
		for (int j=-50; j<=50; j++) {
			if (i==1 && j>0) continue;
			double ratio = 1.0/i+0.001*j;
			ratios.insert(ratio);
		}
	}

	for (set<double>::iterator iter = ratios.begin(); iter != ratios.end(); iter++) {
		double ul = 0.1*util-hDigraph->linear_factor;
		double tl = 1.0*hDigraph->gran_digraph->lper/(*iter);
		
		double wcet = ul*tl;

		cout << "wcet = " << wcet << endl;
		cout << "lfac=" << exactUtil << "\tc_rbf = " << hDigraph->c_rbf << "\tc_ibf = " << hDigraph->c_ibf << endl;
		
		// exact
		double rt = 0;
		for (int t=1; t<500000; t++) {
			rt = ibf_map[t] + wcet;
			if (rt <= t) break;
			t = rt;
			if ( t >= 500000) {
				cerr << "We have arrived at the boundary." << endl;
				exit(EXIT_FAILURE);
			}
		}

		ExactResult.push_back(rt/tl);

		// lubrbf
		rt = (hDigraph->c_rbf+wcet)/(1.0-hDigraph->linear_factor);
		LUBRBFResult.push_back(rt/tl);

		// lubibf
		rt = (hDigraph->c_ibf+wcet)/(1.0-hDigraph->linear_factor);
		LUBIBFResult.push_back(rt/tl);
	}

	delete hDigraph;

}

bool DigraphExperiments::testEffectOfDigraphPeriodWithHarmonicSet(vector<double> &IBFResult, vector<double> &RBFResult, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double>& ratios,  const char *p, int util, double &exactUtil) {
	Digraph* hDigraph;
	FileReader::DotFileReader(hDigraph,1,p);

	// set the utilization
	hDigraph->generate_strongly_connected_components();
	hDigraph->check_strongly_connected();
	hDigraph->calculate_period_gcd();
	hDigraph->calculate_all_gcd();
	hDigraph->calculate_linear_factor();

	double factor = 0.02*util/hDigraph->linear_factor;
	//cout << "Before scaling wcets: factor=" << factor << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << hDigraph->linear_factor << endl;

	// release sccs
	for (vector<Digraph*>::iterator iter = hDigraph->sccs.begin(); iter != hDigraph->sccs.end(); iter++) {
		(*iter)->node_vec.clear();
		(*iter)->edge_vec.clear();
		delete *iter;
		*iter = NULL;
	}
	hDigraph->sccs.clear();

	hDigraph->scale_wcet(factor);
	hDigraph->generate_strongly_connected_components();
	hDigraph->check_strongly_connected();

	hDigraph->calculate_period_gcd();
	hDigraph->calculate_all_gcd();
	hDigraph->calculate_linear_factor();
	hDigraph->calculate_linear_upper_bounds();

	if (hDigraph->strongly_connected) {
		hDigraph->unit_digraph = new UnitDigraph(hDigraph);
		hDigraph->unit_digraph->prepare3(false);
	}
	else {
		cout << "No strongly connected." << endl;
		delete hDigraph;
		return false;
	}

	//cout << "After scaling wcets: factor=" << factor << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << hDigraph->linear_factor << endl;
	
	exactUtil = hDigraph->linear_factor;

	// calculate the periodicity parameters
	//int tf_max = hDigraph->unit_digraph->ldef * hDigraph->pGCD + hDigraph->unit_digraph->lper*hDigraph->pGCD*2000;
	long long int tf_max = (long long int) hDigraph->unit_digraph->lper*hDigraph->pGCD*2000;
	/*
	cout << "tf_max = " << tf_max << endl;
	if (tf_max <= 0) {
		cout << "Too large tf_max!" << endl;
		return false;
	}
	*/
	if (hDigraph->strongly_connected) {
		hDigraph->calculate_rbf_with_periodicity_DP(tf_max,hDigraph->rbf_map3_DP);
		hDigraph->calculate_ibf_with_periodicity_DP(tf_max,hDigraph->ibf_map3_DP);
	}
	else {
		hDigraph->calculate_rbf_without_periodicity_DP(tf_max,hDigraph->rbf_map2_DP);
		hDigraph->calculate_ibf_without_periodicity_DP(tf_max,hDigraph->ibf_map2_DP);
	}

	//cout << "LinearPeriod=" << hDigraph->unit_digraph->lper * hDigraph->pGCD << "\tLinearDefect=" << hDigraph->unit_digraph->ldef*hDigraph->pGCD + 1 << endl;
	//cout << "Crbf = " << hDigraph->c_rbf << "\t Cibf = " << hDigraph->c_ibf << endl;

	for (int i=0; i<1000; i++) {
		double ratio = 1.0-0.001*i;
		ratios.insert(ratio);
	}

	for (set<double>::iterator iter = ratios.begin(); iter != ratios.end(); iter++) {
		double ul = 0.1*util-hDigraph->linear_factor;
		//double tl = 1.0*hDigraph->unit_digraph->lper*hDigraph->pGCD/(*iter) + hDigraph->unit_digraph->ldef * hDigraph->pGCD+1;
		double tl = 1.0*hDigraph->unit_digraph->lper*hDigraph->pGCD/(*iter);
		double wcet = ul*tl;

		//cout << "wcet = " << wcet << endl;
		//cout << "lfac=" << exactUtil << "\tc_rbf = " << hDigraph->c_rbf << "\tc_ibf = " << hDigraph->c_ibf << endl;
		
		// ibf, arbitrary deadlines
		double rt = 0;
		int index = 1;
		while (true) {
			double temp = 0;
			for (long long int t=1; t<tf_max; t++) {
				temp = hDigraph->ibf2(t,hDigraph->ibf_map3_DP) + index * wcet;
				if (temp <= t) break;
				t = temp;
				if ( t >= tf_max) {
					delete hDigraph;
					cout << "We have arrived at the boundary." << endl;
					return false;
				}
			}
			temp = temp - (index-1) * tl;
			if (temp > rt && index > 1) cout << "Index=" << index <<"\trt_ibf=" << temp << endl;
			if (temp <= rt) break;
			rt = max(rt, temp);
			if (rt <= tl) break;
			index ++;
		}

		IBFResult.push_back(rt/tl);
		//cout << "rt_ibf=" << rt << "\tindex=" << index << "\t";

		// rbf, arbitrary deadlines
		rt = 0;
		index = 1;
		while (true) {
			double temp = 0;
			for (long long int t=1; t<tf_max; t++) {
				temp = hDigraph->rbf2(t,hDigraph->rbf_map3_DP) + index * wcet;
				if (temp <= t) break;
				t = temp;
				if ( t >= tf_max) {
					delete hDigraph;
					cout << "We have arrived at the boundary." << endl;
					return false;
				}
			}
			temp = temp - (index-1) * tl;
			if (temp > rt && index > 1)  cout << "Index=" << index <<"\trt_rbf=" << temp << endl;
			if (temp <= rt) break;
			rt = max(rt, temp);
			if (rt <= tl) break;
			index++;
		}

		RBFResult.push_back(rt/tl);
		//cout << "rt_rbf=" << rt << "\tindex=" << index << "\t";

		// lubrbf
		rt = (hDigraph->c_rbf+wcet)/(1.0-hDigraph->linear_factor);
		LUBRBFResult.push_back(rt/tl);
		//cout << rt << "\t";

		// lubibf
		rt = (hDigraph->c_ibf+wcet)/(1.0-hDigraph->linear_factor);
		LUBIBFResult.push_back(rt/tl);
		//cout << rt << endl;
	}

	delete hDigraph;
	return true;
}

bool DigraphExperiments::testEffectOfDigraphPeriodWithHarmonicSetUtil(vector<double> &IBFResult, vector<double> &RBFResult, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double> ratios,  const char *p, int util, double &exactUtil) {
	Digraph* hDigraph;
	FileReader::DotFileReader(hDigraph,1,p);

	// set the utilization
	hDigraph->generate_strongly_connected_components();
	hDigraph->check_strongly_connected();
	hDigraph->calculate_period_gcd();
	hDigraph->calculate_all_gcd();
	hDigraph->calculate_linear_factor();
	hDigraph->calculate_linear_upper_bounds();

	if (hDigraph->strongly_connected) {
		hDigraph->unit_digraph = new UnitDigraph(hDigraph);
		hDigraph->unit_digraph->prepare3(false);
	}
	else {
		cout << "No strongly connected." << endl;
		delete hDigraph;
		return false;
	}

	//cout << "After scaling wcets: factor=" << factor << "\tExpectedUtil=" << 0.1*util << "\tExactUtil=" << hDigraph->linear_factor << endl;
	
	exactUtil = hDigraph->linear_factor;

	// calculate the periodicity parameters
	//int tf_max = hDigraph->unit_digraph->ldef * hDigraph->pGCD + hDigraph->unit_digraph->lper*hDigraph->pGCD*2000;
	long long int tf_max = (long long int) hDigraph->unit_digraph->lper*hDigraph->pGCD*2000;
	/*
	cout << "tf_max = " << tf_max << endl;
	if (tf_max <= 0) {
		cout << "Too large tf_max!" << endl;
		return false;
	}
	*/
	if (hDigraph->strongly_connected) {
		hDigraph->calculate_rbf_with_periodicity_DP(tf_max,hDigraph->rbf_map3_DP);
		hDigraph->calculate_ibf_with_periodicity_DP(tf_max,hDigraph->ibf_map3_DP);
	}
	else {
		hDigraph->calculate_rbf_without_periodicity_DP(tf_max,hDigraph->rbf_map2_DP);
		hDigraph->calculate_ibf_without_periodicity_DP(tf_max,hDigraph->ibf_map2_DP);
	}

	//cout << "LinearPeriod=" << hDigraph->unit_digraph->lper * hDigraph->pGCD << "\tLinearDefect=" << hDigraph->unit_digraph->ldef*hDigraph->pGCD + 1 << endl;
	//cout << "Crbf = " << hDigraph->c_rbf << "\t Cibf = " << hDigraph->c_ibf << endl;

	for (set<double>::iterator iter = ratios.begin(); iter != ratios.end(); iter++) {
		double ul = 0.01*util-hDigraph->linear_factor;
		//double tl = 1.0*hDigraph->unit_digraph->lper*hDigraph->pGCD/(*iter) + hDigraph->unit_digraph->ldef * hDigraph->pGCD+1;
		double tl = 1.0*hDigraph->unit_digraph->lper*hDigraph->pGCD/(*iter);
		double wcet = ul*tl;

		//cout << "wcet = " << wcet << endl;
		//cout << "lfac=" << exactUtil << "\tc_rbf = " << hDigraph->c_rbf << "\tc_ibf = " << hDigraph->c_ibf << endl;
		
		// ibf, arbitrary deadlines
		double rt = 0;
		int index = 1;
		while (true) {
			double temp = 0;
			for (long long int t=1; t<tf_max; t++) {
				temp = hDigraph->ibf2(t,hDigraph->ibf_map3_DP) + index * wcet;
				if (temp <= t) break;
				t = temp;
				if ( t >= tf_max) {
					delete hDigraph;
					cout << "We have arrived at the boundary." << endl;
					return false;
				}
			}
			temp = temp - (index-1) * tl;
			if (temp > rt && index > 1) cout << "Index=" << index <<"\trt_ibf=" << temp << endl;
			if (temp <= rt) break;
			rt = max(rt, temp);
			if (rt <= tl) break;
			index ++;
		}

		IBFResult.push_back(rt/tl);
		//cout << "rt_ibf=" << rt << "\tindex=" << index << "\t";

		// rbf, arbitrary deadlines
		rt = 0;
		index = 1;
		while (true) {
			double temp = 0;
			for (long long int t=1; t<tf_max; t++) {
				temp = hDigraph->rbf2(t,hDigraph->rbf_map3_DP) + index * wcet;
				if (temp <= t) break;
				t = temp;
				if ( t >= tf_max) {
					delete hDigraph;
					cout << "We have arrived at the boundary." << endl;
					return false;
				}
			}
			temp = temp - (index-1) * tl;
			if (temp > rt && index > 1)  cout << "Index=" << index <<"\trt_rbf=" << temp << endl;
			if (temp <= rt) break;
			rt = max(rt, temp);
			if (rt <= tl) break;
			index++;
		}

		RBFResult.push_back(rt/tl);
		//cout << "rt_rbf=" << rt << "\tindex=" << index << "\t";

		// lubrbf
		rt = (hDigraph->c_rbf+wcet)/(1.0-hDigraph->linear_factor);
		LUBRBFResult.push_back(rt/tl);
		//cout << rt << "\t";

		// lubibf
		rt = (hDigraph->c_ibf+wcet)/(1.0-hDigraph->linear_factor);
		LUBIBFResult.push_back(rt/tl);
		//cout << rt << endl;
	}

	delete hDigraph;
	return true;
}

void DigraphExperiments::generateDigraphsForTestingEffectOfTaskParameters(string directory, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int maxRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int exception = 0;
	
	for (int tUtil = minUtil; tUtil <= maxUtil; tUtil += stepUtil) {
		for (int num = minNum; num <= maxNum; num += stepNum) {
			for (int run=0; run<maxRun; run++) {
				
				Digraph** digraphs = new Digraph*[num];
				double* util = Utility::uniformly_distributed(num,0.01*tUtil);
				double totalUtil = 0;
				int times2 = 0;

				try {
					while (1) {
						totalUtil = 0;
						for (int i=0; i<num; i++) {
							Digraph* digraph;
							int times = 0;

							while(1) {
								digraph = RandomGenerator::generate_DIGRAPH_V(i,util[i]);
								digraph->generate_strongly_connected_components();
								digraph->check_strongly_connected();

								if (digraph->strongly_connected) break;
								cout << "times=" << times++ << endl;

								delete digraph;
							}

							digraph->calculate_period_gcd();
							digraph->calculate_all_gcd();
							digraph->calculate_linear_factor();

							totalUtil += digraph->linear_factor;

							/*
							// calculate the linear period
							digraph->gran_digraph = new GranDigraph(digraph);
							digraph->gran_digraph->generate_gran_digraph();
							digraph->gran_digraph->n_size = digraph->gran_digraph->node_vec.size();
							digraph->gran_digraph->generate_exec_request_matrix();
							digraph->gran_digraph->calculate_linear_period(false);

							cout << "lper=" << digraph->gran_digraph->lper << endl;
							*/

							digraphs[i] = digraph;
						}

						cout << times2++ << "=>" << "expected util = " << 0.01*tUtil << "\ttotalUtil = " << totalUtil << endl;

						if (abs(totalUtil-0.01*tUtil) <= 0.01)
							break;

						for (int i=0; i<num; i++)
							delete digraphs[i];
					}

					string name = directory + "\\digraphNum" + Utility::int_to_string(num) + "Util" + Utility::int_to_string(tUtil) +"Run" + Utility::int_to_string(run) + ".dot";
					const char *p = name.c_str();
					cout << name <<endl;
					cout << "expected util = " << 0.01*tUtil << "\ttotalUtil = " << totalUtil << endl;
					FileWriter::DotFileWriter(digraphs, num, p);

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;

					delete[] util;

				}  catch (bad_alloc) {
					exception++;
					cout << "bad_alloc exception. " << exception << endl;
				
					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;

					delete[] util;
				}
			}
		}
	}
}

void DigraphExperiments::testEffectOfTaskParameters(string directory, string file, bool isExact, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun, int choice) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}
	
	int exception = 0;

	map<int,map<int,double>> resultsIBF;
	map<int,map<int,double>> resultsRBF;
	map<int,map<int,double>> resultsLUBRBF;
	map<int,map<int,double>> resultsLUBIBF;
	map<int,map<int,int>> resultsNonSucc;
	
	for (int tUtil = minUtil; tUtil <= maxUtil; tUtil += stepUtil) {
		map<int,double> rIBF;
		map<int,double> rRBF;
		map<int,double> rLUBRBF;
		map<int,double> rLUBIBF;
		map<int,int> rNonSucc;

		for (int num = minNum; num <= maxNum; num += stepNum) {
			int numNonSucc = 0;
			for (int run=startRun; run<=endRun; run++) {
				Digraph** digraphs; 

				string name = directory + "\\digraphNum" + Utility::int_to_string(num) + "Util" + Utility::int_to_string(tUtil) +"Run" + Utility::int_to_string(run) + ".dot";
				const char *p = name.c_str();
				cout << name <<endl;

				int totalNum;
				FileReader::DotFileReader(digraphs, totalNum,1, p);

				if (totalNum != num) {
					cerr << "Error in the task num!" << "\texpectedNum=" <<num << "\trealNum=" <<totalNum << endl;
					exit(EXIT_FAILURE);
				}

				bool success = true;

				try {
					int tf = 5000000;
					int maxPeriod = 1000000;
					// prepare
					for (int i=0; i<num; i++) {
						Digraph* digraph = digraphs[i];
						digraph->calculate_period_gcd();
						digraph->calculate_all_gcd();

						digraph->generate_strongly_connected_components();
						digraph->check_strongly_connected();
						digraph->calculate_linear_factor();

						if (digraph->linear_factor <= 0) {
							success = false;
							goto END;
						}

						digraph->calculate_linear_upper_bounds();

						/*
						// calculate the linear period
						cout << "Calculating the linear period ..." <<endl;
						digraph->unit_digraph = new UnitDigraph(digraph);
						digraph->unit_digraph->generate_unit_digraph();
						digraph->unit_digraph->n_size = digraph->unit_digraph->node_vec.size();
						digraph->unit_digraph->generate_exec_request_matrix();
						digraph->unit_digraph->calculate_linear_period(false);

						maxPeriod = max(maxPeriod, digraph->unit_digraph->lper);
						*/
						

						// calculate rbf and ibf
						cout << "Calculating RBF ..." <<endl;
						digraph->calculate_rbf_without_periodicity_fast(tf,digraph->rbf_map2_fast);
						cout << "Calculating IBF ..." <<endl;
						digraph->calculate_ibf_without_periodicity_fast(tf,digraph->ibf_map2_fast);
					}

					cout << "the maximum period = " << maxPeriod << endl;

					double wcet = 0.1*maxPeriod;
					
					// calculate the response times based on ibf
					double rt_ibf = 0;
					int index = 1;
					while (true) {
						double temp = 0;
						for (int t=1; t<=tf; t++) {
							temp = index*wcet;
							for (int i=0; i<num; i++) 
								temp += digraphs[i]->ibf2_fast(t);
							if (temp <= t) break;
							t = temp;
							if ( t >= tf) {
								cout << "We have arrived at the boundary." << endl;
								success = false;
								goto END;
							}
						}
						temp = temp - (index-1)*maxPeriod;
						if (temp > rt_ibf && index > 1) cout << "Index=" << index <<"\trt_ibf=" << temp << endl;
						if (temp <= rt_ibf) break;
						rt_ibf = max(rt_ibf, temp);
						if (rt_ibf <= maxPeriod) break;
						index ++;
					}

					// calculate the response times based on rbf
					double rt_rbf = 0;
					index = 1;
					while (true) {
						double temp = 0;
						for (int t=1; t<=tf; t++) {
							temp = index*wcet;
							for (int i=0; i<num; i++) 
								temp += digraphs[i]->rbf2_fast(t);
							if (temp <= t) break;
							t = temp;
							if ( t >= tf) {
								cout << "We have arrived at the boundary." << endl;
								success = false;
								goto END;
							}
						}
						temp = temp - (index-1)*maxPeriod;
						if (temp > rt_rbf && index > 1) cout << "Index=" << index <<"\trt_rbf=" << temp << endl;
						if (temp <= rt_rbf) break;
						rt_rbf = max(rt_rbf, temp);
						if (rt_rbf <= maxPeriod) break;
						index ++;
					}

					// calculate the upper bound response times
					double sum_c_rbf = 0;
					double sum_c_ibf = 0;
					double sum_util = 0;

					for (int i=0; i<num; i++) {
						sum_c_rbf += digraphs[i]->c_rbf;
						sum_c_ibf += digraphs[i]->c_ibf;
						sum_util += digraphs[i]->linear_factor;
					}

					double rt_lubrbf = (sum_c_rbf + wcet)/(1.0-sum_util);
					double rt_lubibf = (sum_c_ibf + wcet)/(1.0-sum_util);

					if (rIBF.find(num) == rIBF.end()) {
						rIBF[num] = rt_ibf/maxPeriod;
						rRBF[num] = rt_rbf/maxPeriod;;
						rLUBIBF[num] = rt_lubibf/maxPeriod;
						rLUBRBF[num] = rt_lubrbf/maxPeriod;
					}
					else {
						rIBF[num] += rt_ibf/maxPeriod;
						rRBF[num] += rt_rbf/maxPeriod;;
						rLUBIBF[num] += rt_lubibf/maxPeriod;
						rLUBRBF[num] += rt_lubrbf/maxPeriod;
					}

					cout << "WCET=" << wcet << ":\tRBF = " << rt_rbf << "\tIBF = " << rt_ibf << "\tLUBIBF = " << rt_lubibf << "\tLUBRBF = " << rt_lubrbf << endl;
					fout << "WCET=" << wcet << ":\tRBF = " << rt_rbf << "\tIBF = " << rt_ibf << "\tLUBIBF = " << rt_lubibf << "\tLUBRBF = " << rt_lubrbf << endl;

END:
					if (!success) numNonSucc++;

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;


				} catch (bad_alloc) {
					exception++;
					cout << "bad_alloc exception. " << exception << endl;

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;
				}

				if (success) {
					fout << "====================Util=" << tUtil << ",Num=" << num << ", run=" << run+1 << ", nonSucc=" << numNonSucc << "=======================" << endl; 
					outputResults(fout, rIBF, "IBF"+Utility::int_to_string(tUtil));
					outputResults(fout, rRBF, "RBF"+Utility::int_to_string(tUtil));
					outputResults(fout, rLUBRBF, "LUBRBF"+Utility::int_to_string(tUtil));
					outputResults(fout, rLUBIBF, "LUBIBF"+Utility::int_to_string(tUtil));
				}
			}

			rNonSucc[num] = numNonSucc;
		}

		resultsIBF[tUtil] = rIBF;
		resultsRBF[tUtil]= rRBF;
		resultsLUBRBF[tUtil] = rLUBRBF;
		resultsLUBIBF[tUtil] = rLUBIBF;
		resultsNonSucc[tUtil] = rNonSucc;
	}

	fout << "===============ALL=================" << endl;
	outputResults(fout, resultsIBF, "IBF");
	outputResults(fout, resultsRBF, "RBF");
	outputResults(fout, resultsLUBRBF, "LUBRBF");
	outputResults(fout, resultsLUBIBF, "LUBIBF");
	outputResults(fout, resultsNonSucc, "numNonSucc");

	fout.close();
}

void DigraphExperiments::generateDigraphsForTestingEffectOfTaskParameters2(string directory, int maxNodeNum, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun) {
	// mkdir directory
	if ( _mkdir(directory.c_str()) == 0) {
		cout<<"Directory: "<<directory<<" was successfully created."<<endl;
	}

	int exception = 0;
	
	for (int tUtil = minUtil; tUtil <= maxUtil; tUtil += stepUtil) {
		for (int num = minNum; num <= maxNum; num += stepNum) {
			for (int run=startRun; run<=endRun; run++) {
				Digraph** digraphs;

				try {
					digraphs = RandomGenerator::generate_digraph_system_sccs(num,maxNodeNum,0.01*tUtil,false);

					double realUtil = 0;
					for (int i=0; i<num; i++) {
						Digraph* digraph = digraphs[i];
						realUtil += digraph->linear_factor;
					}

					cout << "ExpectedUtil=" << 0.01*tUtil << "\tRealUtil=" << realUtil << endl;

					string name = directory + "\\digraphNum" + Utility::int_to_string(num) + "Util" + Utility::int_to_string(tUtil) +"Run" + Utility::int_to_string(run) + ".dot";
					const char *p = name.c_str();
					cout << name <<endl;
					FileWriter::DotFileWriter(digraphs, num, p);

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;

				}  catch (bad_alloc) {
					exception++;
					cout << "bad_alloc exception. " << exception << endl;
				
					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;
				}
			}
		}
	}
}

void DigraphExperiments::testEffectOfTaskParameters2(string directory, string file, double factor, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun, int choice) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}
	
	int exception = 0;

	map<int,map<int,double>> resultsIBF;
	map<int,map<int,double>> resultsRBF;
	map<int,map<int,double>> resultsLUBRBF;
	map<int,map<int,double>> resultsLUBIBF;
	map<int,map<int,int>> resultsNonSucc;
	
	for (int tUtil = minUtil; tUtil <= maxUtil; tUtil += stepUtil) {
		map<int,double> rIBF;
		map<int,double> rRBF;
		map<int,double> rLUBRBF;
		map<int,double> rLUBIBF;
		map<int,int> rNonSucc;

		for (int num = minNum; num <= maxNum; num += stepNum) {
			int numNonSucc = 0;

			for (int run=startRun; run<=endRun; run++) {
				Digraph** digraphs; 

				string name = directory + "\\digraphNum" + Utility::int_to_string(num) + "Util" + Utility::int_to_string(tUtil) +"Run" + Utility::int_to_string(run) + ".dot";
				const char *p = name.c_str();
				cout << name <<endl;

				int totalNum;
				FileReader::DotFileReader(digraphs, totalNum,1, p);

				if (totalNum != num) {
					cerr << "Error in the task num!" << "\texpectedNum=" <<num << "\trealNum=" <<totalNum << endl;
					exit(EXIT_FAILURE);
				}

				bool success = true;

				try {
					long long int maxPeriod = INT_MIN;
					// prepare
					for (int i=0; i<num; i++) {
						Digraph* digraph = digraphs[i];
						digraph->calculate_period_gcd();
						digraph->calculate_all_gcd();

						digraph->generate_strongly_connected_components();
						digraph->check_strongly_connected();
						digraph->calculate_linear_factor();

						if (digraph->linear_factor <= 0  || !digraph->strongly_connected) {
							success = false;
							goto END;
						}

						digraph->calculate_linear_upper_bounds();

						digraph->unit_digraph = new UnitDigraph(digraph);
						digraph->unit_digraph->prepare3(false);

						maxPeriod = max(maxPeriod, digraph->unit_digraph->lper*digraph->pGCD);
					}

					//maxPeriod = (double)maxPeriod/(0.01*factor)+linearDefect;
					maxPeriod = (double)maxPeriod/(0.01*factor);

					//double wcet = 0.01*4*tUtil*maxPeriod;
					double wcet = 0.1*maxPeriod;

					// calculate the upper bound response times
					double sum_c_rbf = 0;
					double sum_c_ibf = 0;
					double sum_util = 0;

					for (int i=0; i<num; i++) {
						sum_c_rbf += digraphs[i]->c_rbf;
						sum_c_ibf += digraphs[i]->c_ibf;
						sum_util += digraphs[i]->linear_factor;
					}

					double rt_lubrbf = (sum_c_rbf + wcet)/(1.0-sum_util);
					double rt_lubibf = (sum_c_ibf + wcet)/(1.0-sum_util);

					long long int tf = rt_lubrbf;
					//if (tf <=0) tf = INT_MAX;
					cout << "the maximum period = " << maxPeriod << "\ttf=" << tf << endl;

					int linearDefect = INT_MIN;
					for (int i=0; i<num; i++) {
						Digraph* digraph = digraphs[i];
						// calculate rbf and ibf
						cout << "Calculating RBF ..." <<endl;
						digraph->calculate_rbf_with_periodicity_DP(tf,digraph->rbf_map3_DP);
						cout << "Calculating IBF ..." <<endl;
						digraph->calculate_ibf_with_periodicity_DP(tf,digraph->ibf_map3_DP);
						linearDefect = max(linearDefect, digraph->unit_digraph->ldef * digraph->pGCD + 1);
					}
					
					// calculate the response times based on ibf
					double rt_ibf = 0;
					int index = 1;
					while(true) {
						double temp =0;
						for (long long int t=rt_ibf+1; t<=tf; t++) {
							temp = wcet*index;
							for (int i=0; i<num; i++) 
								temp += digraphs[i]->ibf2(t,digraphs[i]->ibf_map3_DP);
							if (temp <= t) break;
							t = temp;
							if ( t >= tf) {
								cout << "We have arrived at the boundary." << endl;
								success = false;
								goto END;
							}
						}
						temp = temp - (index-1)*maxPeriod;
						if (temp > rt_ibf && index > 1) 
							cout << "Index=" << index <<"\trt_ibf=" << temp << endl;
						if (temp <= rt_ibf) break;
						rt_ibf = max(rt_ibf, temp);
						if (rt_ibf <= maxPeriod) break;
						index ++;
					}

					// calculate the response times based on rbf
					double rt_rbf = 0;
					index = 1;
					while(true) {
						double temp =0;
						for (long long int t=rt_rbf+1; t<=tf; t++) {
							temp = wcet*index;
							for (int i=0; i<num; i++) 
								temp += digraphs[i]->rbf2(t,digraphs[i]->rbf_map3_DP);
							if (temp <= t) break;
							t = temp;
							if ( t >= tf) {
								cout << "We have arrived at the boundary." << endl;
								success = false;
								goto END;
							}
						}
						temp = temp - (index-1)*maxPeriod;
						if (temp > rt_rbf && index > 1) 
							cout << "Index=" << index <<"\trt_rbf=" << temp << endl;
						if (temp <= rt_rbf) break;
						rt_rbf = max(rt_rbf, temp);
						if (rt_rbf <= maxPeriod) break;
						index ++;
					}

					if (rIBF.find(num) == rIBF.end()) {
						rIBF[num] = rt_ibf/maxPeriod;
						rRBF[num] = rt_rbf/maxPeriod;;
						rLUBIBF[num] = rt_lubibf/maxPeriod;
						rLUBRBF[num] = rt_lubrbf/maxPeriod;
					}
					else {
						rIBF[num] += rt_ibf/maxPeriod;
						rRBF[num] += rt_rbf/maxPeriod;;
						rLUBIBF[num] += rt_lubibf/maxPeriod;
						rLUBRBF[num] += rt_lubrbf/maxPeriod;
					}

					cout << "WCET=" << wcet << ":\tRBF = " << rt_rbf << "\tIBF = " << rt_ibf << "\tLUBIBF = " << rt_lubibf << "\tLUBRBF = " << rt_lubrbf << endl;
					fout << "WCET=" << wcet << ":\tRBF = " << rt_rbf << "\tIBF = " << rt_ibf << "\tLUBIBF = " << rt_lubibf << "\tLUBRBF = " << rt_lubrbf << endl;

END:
					if (!success) numNonSucc++;

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;


				} catch (bad_alloc) {
					exception++;
					cerr << "bad_alloc exception. " << exception << endl;
					exit(EXIT_FAILURE);

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;
				}

				if (success) {
					fout << "====================Util=" << tUtil << ",Num=" << num << ", run=" << run+1 << ", nonSucc=" << numNonSucc << "=======================" << endl; 
					outputResults(fout, rIBF, "IBF"+Utility::int_to_string(tUtil));
					outputResults(fout, rRBF, "RBF"+Utility::int_to_string(tUtil));
					outputResults(fout, rLUBRBF, "LUBRBF"+Utility::int_to_string(tUtil));
					outputResults(fout, rLUBIBF, "LUBIBF"+Utility::int_to_string(tUtil));
				}
			}

			rNonSucc[num] = numNonSucc;
		}

		resultsIBF[tUtil] = rIBF;
		resultsRBF[tUtil]= rRBF;
		resultsLUBRBF[tUtil] = rLUBRBF;
		resultsLUBIBF[tUtil] = rLUBIBF;
		resultsNonSucc[tUtil] = rNonSucc;
	}

	fout << "===============ALL=================" << endl;
	outputResults(fout, resultsIBF, "IBF");
	outputResults(fout, resultsRBF, "RBF");
	outputResults(fout, resultsLUBRBF, "LUBRBF");
	outputResults(fout, resultsLUBIBF, "LUBIBF");
	outputResults(fout, resultsNonSucc, "numNonSucc");

	fout.close();
}

void DigraphExperiments::testEffectOfTaskParameters2Util(string directory, string file, double factor, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun, int choice) {
	ofstream fout(file, ios::out | ios::trunc);
	if (!fout.is_open()) {
		cerr << "Can't open "<<file<<" file for output" <<endl;
		exit(EXIT_FAILURE);
	}

	int exception = 0;

	map<int,map<int,double>> resultsIBF;
	map<int,map<int,double>> resultsRBF;
	map<int,map<int,double>> resultsLUBRBF;
	map<int,map<int,double>> resultsLUBIBF;
	map<int,map<int,int>> resultsNonSucc;

	for (int num = minNum; num <= maxNum; num += stepNum) {
		map<int,double> rIBF;
		map<int,double> rRBF;
		map<int,double> rLUBRBF;
		map<int,double> rLUBIBF;
		map<int,int> rNonSucc;

		for (int tUtil = minUtil; tUtil <= maxUtil; tUtil += stepUtil) {
			int numNonSucc = 0;

			for (int run=startRun; run<=endRun; run++) {
				Digraph** digraphs; 

				string name = directory + "\\digraphNum" + Utility::int_to_string(num) + "Util" + Utility::int_to_string(tUtil) +"Run" + Utility::int_to_string(run) + ".dot";
				const char *p = name.c_str();
				cout << name <<endl;

				int totalNum;
				FileReader::DotFileReader(digraphs, totalNum,1, p);

				if (totalNum != num) {
					cerr << "Error in the task num!" << "\texpectedNum=" <<num << "\trealNum=" <<totalNum << endl;
					exit(EXIT_FAILURE);
				}

				bool success = true;

				try {
					long long int maxPeriod = INT_MIN;
					// prepare
					for (int i=0; i<num; i++) {
						Digraph* digraph = digraphs[i];
						digraph->calculate_period_gcd();
						digraph->calculate_all_gcd();

						digraph->generate_strongly_connected_components();
						digraph->check_strongly_connected();
						digraph->calculate_linear_factor();

						if (digraph->linear_factor <= 0  || !digraph->strongly_connected) {
							success = false;
							goto END;
						}

						digraph->calculate_linear_upper_bounds();

						digraph->unit_digraph = new UnitDigraph(digraph);
						digraph->unit_digraph->prepare3(false);

						maxPeriod = max(maxPeriod, digraph->unit_digraph->lper*digraph->pGCD);
					}

					//maxPeriod = (double)maxPeriod/(0.01*factor)+linearDefect;
					maxPeriod = (double)maxPeriod/(0.01*factor);

					//double wcet = 0.01*4*tUtil*maxPeriod;
					double wcet = 0.1*maxPeriod;

					// calculate the upper bound response times
					double sum_c_rbf = 0;
					double sum_c_ibf = 0;
					double sum_util = 0;

					for (int i=0; i<num; i++) {
						sum_c_rbf += digraphs[i]->c_rbf;
						sum_c_ibf += digraphs[i]->c_ibf;
						sum_util += digraphs[i]->linear_factor;
					}

					double rt_lubrbf = (sum_c_rbf + wcet)/(1.0-sum_util);
					double rt_lubibf = (sum_c_ibf + wcet)/(1.0-sum_util);

					long long int tf = rt_lubrbf;
					//if (tf <=0) tf = INT_MAX;
					cout << "the maximum period = " << maxPeriod << "\ttf=" << tf << endl;

					int linearDefect = INT_MIN;
					for (int i=0; i<num; i++) {
						Digraph* digraph = digraphs[i];
						// calculate rbf and ibf
						cout << "Calculating RBF ..." <<endl;
						digraph->calculate_rbf_with_periodicity_DP(tf,digraph->rbf_map3_DP);
						cout << "Calculating IBF ..." <<endl;
						digraph->calculate_ibf_with_periodicity_DP(tf,digraph->ibf_map3_DP);
						linearDefect = max(linearDefect, digraph->unit_digraph->ldef * digraph->pGCD + 1);
					}

					// calculate the response times based on ibf
					double rt_ibf = 0;
					int index = 1;
					while(true) {
						double temp =0;
						for (long long int t=rt_ibf+1; t<=tf; t++) {
							temp = wcet*index;
							for (int i=0; i<num; i++) 
								temp += digraphs[i]->ibf2(t,digraphs[i]->ibf_map3_DP);
							if (temp <= t) break;
							t = temp;
							if ( t >= tf) {
								cout << "We have arrived at the boundary." << endl;
								success = false;
								goto END;
							}
						}
						temp = temp - (index-1)*maxPeriod;
						if (temp > rt_ibf && index > 1) 
							cout << "Index=" << index <<"\trt_ibf=" << temp << endl;
						if (temp <= rt_ibf) break;
						rt_ibf = max(rt_ibf, temp);
						if (rt_ibf <= maxPeriod) break;
						index ++;
					}

					// calculate the response times based on rbf
					double rt_rbf = 0;
					index = 1;
					while(true) {
						double temp =0;
						for (long long int t=rt_rbf+1; t<=tf; t++) {
							temp = wcet*index;
							for (int i=0; i<num; i++) 
								temp += digraphs[i]->rbf2(t,digraphs[i]->rbf_map3_DP);
							if (temp <= t) break;
							t = temp;
							if ( t >= tf) {
								cout << "We have arrived at the boundary." << endl;
								success = false;
								goto END;
							}
						}
						temp = temp - (index-1)*maxPeriod;
						if (temp > rt_rbf && index > 1) 
							cout << "Index=" << index <<"\trt_rbf=" << temp << endl;
						if (temp <= rt_rbf) break;
						rt_rbf = max(rt_rbf, temp);
						if (rt_rbf <= maxPeriod) break;
						index ++;
					}

					if (rIBF.find(tUtil) == rIBF.end()) {
						rIBF[tUtil] = rt_ibf/maxPeriod;
						rRBF[tUtil] = rt_rbf/maxPeriod;;
						rLUBIBF[tUtil] = rt_lubibf/maxPeriod;
						rLUBRBF[tUtil] = rt_lubrbf/maxPeriod;
					}
					else {
						rIBF[tUtil] += rt_ibf/maxPeriod;
						rRBF[tUtil] += rt_rbf/maxPeriod;;
						rLUBIBF[tUtil] += rt_lubibf/maxPeriod;
						rLUBRBF[tUtil] += rt_lubrbf/maxPeriod;
					}

					cout << "WCET=" << wcet << ":\tRBF = " << rt_rbf << "\tIBF = " << rt_ibf << "\tLUBIBF = " << rt_lubibf << "\tLUBRBF = " << rt_lubrbf << endl;
					fout << "WCET=" << wcet << ":\tRBF = " << rt_rbf << "\tIBF = " << rt_ibf << "\tLUBIBF = " << rt_lubibf << "\tLUBRBF = " << rt_lubrbf << endl;

END:
					if (!success) numNonSucc++;

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;


				} catch (bad_alloc) {
					exception++;
					cerr << "bad_alloc exception. " << exception << endl;
					exit(EXIT_FAILURE);

					for (int i=0; i<num; i++)
						delete digraphs[i];
					delete[] digraphs;
				}

				if (success) {
					fout << "====================Util=" << tUtil << ",Num=" << num << ", run=" << run+1 << ", nonSucc=" << numNonSucc << "=======================" << endl; 
					outputResults(fout, rIBF, "IBF"+Utility::int_to_string(num));
					outputResults(fout, rRBF, "RBF"+Utility::int_to_string(num));
					outputResults(fout, rLUBRBF, "LUBRBF"+Utility::int_to_string(num));
					outputResults(fout, rLUBIBF, "LUBIBF"+Utility::int_to_string(num));
				}
			}

			rNonSucc[tUtil] = numNonSucc;
		}

		resultsIBF[num] = rIBF;
		resultsRBF[num]= rRBF;
		resultsLUBRBF[num] = rLUBRBF;
		resultsLUBIBF[num] = rLUBIBF;
		resultsNonSucc[num] = rNonSucc;
	}

	fout << "===============ALL=================" << endl;
	outputResults(fout, resultsIBF, "IBF");
	outputResults(fout, resultsRBF, "RBF");
	outputResults(fout, resultsLUBRBF, "LUBRBF");
	outputResults(fout, resultsLUBIBF, "LUBIBF");
	outputResults(fout, resultsNonSucc, "numNonSucc");

	fout.close();
}

void DigraphExperiments::outputResults(ofstream& fout, vector<vector<double>> results, string name) {
	
	//fout << "============" << name << "==========" <<endl;

	for (int i=0; i<results.size(); i++) {
		vector<double> result = results.at(i);
		fout << name <<  i << "=[";
		for (int j=0; j<result.size(); j++) {
			fout << result.at(j);
			if (j != result.size()-1)
				fout << ",";
		}
		fout << "];" << endl;
	}
	//fout << "=====================================" << endl;
}

void DigraphExperiments::outputResults2(ofstream& fout, vector<vector<double>> results, string name) {

	//fout << "============" << name << "==========" <<endl;
	int sizeA = results.size();
	int sizeB = results.front().size();

	for (int i=0; i<sizeB; i++) {
		fout << name <<  i << "=[";
		for (int j=0; j<sizeA; j++) {
			fout << results[j][i];
			if (j != sizeA-1)
				fout << ",";
		}
		fout << "];" << endl;
	}
	//fout << "=====================================" << endl;
}

void DigraphExperiments::outputResults(ofstream& fout, map<int,map<int,double>> results, string name) {
	
	//fout << "============" << name << "==========" <<endl;

	for (map<int,map<int,double>>::iterator iter = results.begin(); iter != results.end(); iter ++) {
		outputResults(fout,iter->second,name+Utility::int_to_string(iter->first));
	}

	//fout << "=====================================" << endl;
}

void DigraphExperiments::outputResults(ofstream& fout, map<int,map<int,int>> results, string name) {
	
	//fout << "============" << name << "==========" <<endl;

	for (map<int,map<int,int>>::iterator iter = results.begin(); iter != results.end(); iter ++) {
		outputResults(fout,iter->second,name+Utility::int_to_string(iter->first));
	}

	//fout << "=====================================" << endl;
}

void DigraphExperiments::outputResults(ofstream& fout, map<int,double> results, string name) {
	
	//fout << "============" << name << "==========" <<endl;

	fout << name << "=[";
	for (map<int,double>::iterator iter = results.begin(); iter != results.end(); iter ++) {
		fout << iter->second;
		if (iter != (--results.end()))
			fout << ",";
	}
	fout << "]" << endl;
	//fout << "=====================================" << endl;
}

void DigraphExperiments::outputResults(ofstream& fout, map<int,int> results, string name) {
	
	//fout << "============" << name << "==========" <<endl;

	fout << name << "=[";
	for (map<int,int>::iterator iter = results.begin(); iter != results.end(); iter ++) {
		fout << iter->second;
		if (iter != (--results.end()))
			fout << ",";
	}
	fout << "]" << endl;
	//fout << "=====================================" << endl;
}