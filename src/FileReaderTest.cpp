#include "Definitions.h"

#ifdef GOOGLETEST
#include <gtest/gtest.h>

#include "FileReader.h"
#include "SchedulabilityAnalysis.h"

extern bool output;
/*
TEST(FileReaderTest, Test0)
{
	Stateflow** sfs;
	int num;
	const char* file = "Output\\Test0Run886.dot";
	FileReader::DotFileReader(sfs, num,1000,file);
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		//sf->write_graphviz(cout);
		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		// show rbfs
		if (output) {
			for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
				cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
				Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
			}
		}

		sf->generate_exec_req_matrix();
		// show execution request matrix
		if (output) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();
	}
	SchedulabilityAnalysis::prepare_all_stateflows(sfs,num,0);
	bool isRBF = SchedulabilityAnalysis::rbf_analysis_static_offset(sfs, num,0);
	bool isIBF = SchedulabilityAnalysis::ibf_analysis_static_offset(sfs,num,0);

	cout<<isRBF<<"\t"<<isIBF<<endl;
}
*/
/*
TEST(FileReaderTest, Test1)
{
	Digraph** digraphs;
	int num;
	const char* file = "Input\\00\\run0.dot";
	FileReader::DotFileReader(digraphs, num, 1, file);

	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];

		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
		digraph->calculate_linear_factor2();
	}

	SchedulabilityAnalysis::prepare_all_digraphs(digraphs,num);
	bool isRBF = SchedulabilityAnalysis::rbf_analysis(digraphs,num,0);
	cout<<isRBF<<endl;
}

TEST(FileReaderTest, Test2)
{
	Digraph** digraphs;
	int num;
	const char* file = "Input\\01\\run25.dot";
	FileReader::DotFileReader(digraphs, num, 1, file);

	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];

		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
		digraph->calculate_linear_factor();
	}

	SchedulabilityAnalysis::prepare_all_digraphs(digraphs,num);
	bool isRBF = SchedulabilityAnalysis::rbf_analysis(digraphs,num,0);
	cout<<isRBF<<endl;
}

TEST(FileReaderTest, Test3)
{
	Digraph** digraphs;
	int num;
	const char* file = "Input\\11\\run1.dot";
	FileReader::DotFileReader(digraphs, num, 1, file);

	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];

		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
		digraph->calculate_linear_factor();
	}

	SchedulabilityAnalysis::prepare_all_digraphs(digraphs,num);
	bool isRBF = SchedulabilityAnalysis::rbf_analysis(digraphs,num,0);
	cout<<isRBF<<endl;
}
*/

#endif