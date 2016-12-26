/* \file RandomGenerator.cpp
*  this file implements a random generator for Digraph or Stateflow systemss. 
*  \author Chao Peng
*  
*  Changes
*  ------
*  04-Sept-2015 : initial revision (CP)
*
*/

#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_

//#include <vld.h>

#include "Digraph.h"
#include "Stateflow.h"

using namespace std;

class RandomGenerator {
public:
	/// Graph properties
	static int DIGRAPH_SCALE;
	static int DIGRAPH_PERIOD[3][2]; // {{2,4},{6,12},{5,10}};
	static int DIGRAPH_SCALE_FACTOR[5]; // {1,2,4,5,10};
	static int DIGRAPH_SCALE_FACTOR2[4]; // {1,10,100,1000}

	/// Stateflow properties
	static int STATEFLOW_SCALE;
	static enum PERIOD_CHOICE
	{
		ExactAnalysis0, // {10, 20, 25, 50} and {1,2}
		ExactAnalysis1, // {5, 10, 20, 25, 50} and {1,2,4}
		ExactAnalysis2, // {5, 10, 20, 25, 50} and {1,2}
		ExactAnalysis3, // {10, 20, 25, 50} and {1,2,4}
		ExactAnalysis4, // {10, 20, 25, 50}
		ExactAnalysis5, // {10, 20, 40}
		ExactAnalysis6, // {10, 20, 25, 40, 50, 100}
		ExactAnalysis7, // {10, 20, 40, 80}
		ExactAnalysis8, // {5, 10, 20, 25, 50, 100}
		ExactAnalysis9, // {10, 20, 50, 100}
		ApproximateAnalysis0, // {5, 10, 20, 25, 50, 100}
		ApproximateAnalysis1, // {5, 10, 20, 25, 50, 100} and {1,2,4,5}
		ApproximateAnalysis2, // {5, 10, 20, 25, 50, 100} and {1,2,4}
		ApproximateAnalysis3, // {1, 2, 5, 10}^2
		ApproximateAnalysis4, // base {1,2,5,10} and factor {1,2,5,10}
		ApproximateAnalysis5, // base {1,2,5,10,20,25,50,100} and factor {1,2,5,10}
		ApproximateAnalysis6, // {10,20,25,50}
		ApproximateAnalysis7, // {5,10,20,25,50}
		ApproximateAnalysis8, // base {1,2,5,10,20,25,50,100,200,500} and factor {1,2,5,10}
		ApproximateAnalysis9, // base {1,2,5,10,20,25,50,100,200,250} and factor {1,2,5,10}
		ApproximateAnalysis10, // base {1,2,5,10,20,25,50,100,200,1000} and factor {1,2,5,10}
		ApproximateAnalysis11, // base {1,2,5,10,20,25,50,100,1000} and factor {1,2,5,10}
		ApproximateAnalysis12, // base {1,2,5,10,20,50,100,200,1000} and factor {1,2,5,10}
		ApproximateAnalysis13, // base {1,2,5,10,20,25,50,100,1000,2000} and factor {1,2,5,10}
	};

	static int STATEFLOW_PERIOD[6]; // {5, 10, 20, 25, 50, 100};
	static int STATEFLOW_PERIOD2[5]; // {10, 20, 25, 50, 100}
	static int STATEFLOW_PERIOD3[4]; // {10, 20, 25, 50}
	static int STATEFLOW_PERIOD4[5]; // {5, 10, 20, 25, 50}
	static int STATEFLOW_PERIOD5[4]; // {1, 2, 5, 10}
	static int STATEFLOW_PERIOD6[8]; // {1,2,5,10,20,25,50,100}
	static int STATEFLOW_PERIOD7[3]; // {10, 20, 25}
	static int STATEFLOW_PERIOD8[10]; // {1,2,5,10,20,25,50,100,200,500}
	static int STATEFLOW_PERIOD9[10]; // {1,2,5,10,20,25,50,100,200,250}
	static int STATEFLOW_PERIOD10[6]; // {10, 20, 25, 40, 50, 100}
	static int STATEFLOW_PERIOD11[4]; // {10, 20, 40, 80}
	static int STATEFLOW_PERIOD12[10]; // {1,2,5,10,20,25,50,100,200,1000}
	static int STATEFLOW_PERIOD13[9]; // {1,2,5,10,20,25,50,100,1000}
	static int STATEFLOW_PERIOD14[9]; // {1,2,5,10,20,50,100,200,1000}
	static int STATEFLOW_PERIOD15[5]; // {10,20,40,50,100}
	static int STATEFLOW_PERIOD16[4]; // {10,20,50,100}
	static int STATEFLOW_PERIOD17[10]; // {1,2,5,10,20,50,100,200,1000,10000}
	//static int STATEFLOW_FACTOR[6]; // {1, 2, 4, 5, 10, 20}
	static int STATEFLOW_BASE[8]; // {1100, 1300,2300, 3100, 3700, 3000, 5000, 7000}
	//static int STATEFLOW_BASE[5]; // {1100, 1300,2300, 3100, 3700, 3000, 5000}
	static int STATEFLOW_FACTOR[5]; // {1, 2, 4, 5, 10, 20}
	static int STATEFLOW_FACTOR2[4]; // {1,2,4,5}
	static int STATEFLOW_FACTOR3[3]; // {1,2,4}
	static int STATEFLOW_FACTOR4[2]; // {1,2}

	/// Generate one random digraph with nNode nodes and nEdge edges
	/// which is the same as the paper ECRTS2013
	static Digraph* generate_one_digraph(int index, int scale, int nNode, int nEdge);
	static Digraph* generate_one_digraph(int index, int scale, int nNode);
	static Digraph generate_one_digraph_without_pointer(int index, int scale, int nNode);
	static Digraph* generate_one_digraph(int index, int scale, int nNode, double util);
	/// Generate one random digraph with nNode
	static Digraph* generate_one_digraph2(int index, int scale, int nNode);
	/// generate one random digraph similar to the digraph setting in the RTSS2014 paper
	static Digraph* generate_one_digraph3(int index, int scale);
	/// generate one random digraph with the given utilization
	/// Node number = [5,10]
    /// Outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// 100% strongly connected
    /// Inter-release separation time = [10,100]
    /// Deadline = the minimum inter-release time of all the out-going edges
	static Digraph* generate_one_digraph4(int index, int scale, double util);

	/// generate DIGRAPH I with the given utilization (not exact)
	/// vertices = [5,10]
	/// in-out degreee = 3.8
	/// outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// strongly connected
	/// period = [20,100]
	/// scale = 1
	static Digraph* generate_DIGRAPH_I(int index, double util);

	/// generate DIGRAPH II with the given utilization (not exact)
	/// vertices = [5,10]
	/// in-out degreee = 3.8
	/// outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// strongly connected
	/// period = [100,300]*scale
	/// scale = {10,100,1000}
	static Digraph* generate_DIGRAPH_II(int index, double util);

	/// generate DIGRAPH III
	/// vertices = [1,10]
	/// in-out degreee = 3.8
	/// outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// strongly connected
	/// period = [50,100]
	/// wcet = [1,8]
	static Digraph* generate_DIGRAPH_III(int index);

	/// generate DIGRAPH IV
	/// /// vertices = [5,10]
	/// in-out degreee = 3.8
	/// outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// strongly connected
	/// period = [50,100]
	/// wcet = util*period
	static Digraph* generate_DIGRAPH_IV(int index, double util);

	/// generate DIGRAPH V
	/// /// vertices = [1,10]
	/// in-out degreee = 3.8
	/// outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// strongly connected
	/// period = [300,500]
	/// wcet = util*period
	static Digraph* generate_DIGRAPH_V(int index, double util);

	/// generate DIGRAPH VI
	/// /// vertices = [7,15]
	/// in-out degreee = 3.8
	/// outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// strongly connected
	/// two types of verticies:
	/// Type-1: wcet = [1,6] and period = [20,100]
	/// Type-2: wcet = [7,9] and period = [101,300]
	static Digraph* generate_DIGRAPH_VI(int index);

	/// generate DIGRAPH VII
	/// vertices = [5,7]
	/// in-out degreee = 3.8
	/// outgoing edge number = {40% with one edge, 40% with two edges, 10% with three edges and 10% with four edges}
    /// strongly connected
	/// wcet = [1,8]
	/// Type-1: period = [20,40]
	/// Type-2: and period = [40,60]
	/// Type-3: and period = [60,80]
	static Digraph* generate_DIGRAPH_VII(int index, int periodType);

	/// Generate one random digraph systems with num digraphs and tUtil (total) utilization, 
	/// similar to the RTSJ2015 paper
	static Digraph** generate_digraph_system(int num, int maxNode, double tUtil, bool deadlineProperty);
	static Digraph** generate_digraph_system(int num, int minNode, int maxNode, double tUtil, bool deadlineProperty);
	static Digraph** generate_digraph_system_without_pointer(int num, int minNode, int maxNode, double tUtil, bool deadlineProperty);
	/// Generate one random digraph systems with num digraphs and tUtil (total) utilization, 
	/// similar to the RTSJ2015 paper, but all the digraphs are strongly connected 
	static Digraph** generate_digraph_system_sccs(int num, int maxNode, double tUtil, bool deadlineProperty);
	/// Generate one random digraph for the given utilization
	static Digraph** generate_digraph_system2(int& num, double tUtil);
	/// generate one digraph system with the given utilization
	static Digraph** generate_digraph_system3(int num, double tUtil);
	/// generate one digraph system with the given utilization
	static Digraph** generate_digraph_system4(int& num, double tUtil);
	/// generate one digraph system similar to the RTSS2014 paper with the given utilization
	static Digraph** generate_digraph_system_guan(int& num, double tUtil);

	/// generate one digraph system with two tasks.
	/// one is a digraph and another is a periodic task
	//static Digraph** generate_digraph_system5(double tUtil);

	/// generate one digraph system with 50% DIGRAPH I and 50% DIGRAPH II
	static Digraph** generate_mixed_digraph_system(int num, double tUtil);

	static int calculate_base();
	static int calculate_separation_time(int base);

	/// Generate one random stateflow with nState states and nTran transitions
	/// which is the same as the paper ECRTS2012
	static Stateflow* generate_one_stateflow(int index, int scale, int nState, int nTran);
	/// Generate one random stateflow with nState states and nTran transitions
	/// and make sure the out-degree of any state <= 2
	static Stateflow* generate_one_stateflow2(int index, int scale, int nState, int nTran);

	/// Generate one random stateflow with nState states
	/// Force to be strongly connected
	static Stateflow* generate_one_stateflow3(int index, int scale, int nState);

	static void setup_periods_for_one_stateflow_for_exact_analysis_0(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_1(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_2(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_3(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_4(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_5(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_6(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_7(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_8(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_exact_analysis_9(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_0(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_1(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_2(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_3(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_4(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_5(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_6(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_7(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_8(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_9(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_10(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_11(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_12(Stateflow* sf, int scale);
	static void setup_periods_for_one_stateflow_for_approximate_analysis_13(Stateflow* sf, int scale);

	static void setup_wcets_for_one_stateflow(Stateflow* sf, double util);
	static void setup_wcets_for_one_stateflow(Stateflow* sf, double util, double& diffUtil);
	static void setup_wcets_for_one_stateflow2(Stateflow* sf, double util);
	

	static Stateflow* generate_one_stateflow_with_util(int index, int scale, double util, int nState, int nTran);
	static Stateflow* generate_one_stateflow_with_util2(int index, int scale, double util, int nState, int nTran);
	/// Generate one random stateflow system with num stateflows and tUtil (total) utilization
	static Stateflow** generate_stateflow_system(int num, int maxState, double tUtil);

	/// Generate one random stateflow system based on STATEFLOW_FACTOR and STATEFLOW_BASE
	static Stateflow** generate_stateflow_system2(int num, int maxState, double tUtil);

	/// Generate one random stateflow system (two partitions: 1~num-1 and num)
	static Stateflow** generate_stateflow_system3(int num, int maxState, double tUtil);

	/// Generate one random stateflow system (the last one is special)
	static Stateflow** generate_stateflow_system4(int num, int maxState, double tUtil);
	static Stateflow** generate_stateflow_system5(int num, int maxState, double tUtil);
	/// if period_choice == 1, we use the function setup_periods_for_one_stateflow
	/// else if period_choice == 2, we use the function setup_periods_for_one_stateflow2
	/// else if period_choice == 3, we use the function setup_periods_for_one_stateflow3
	static Stateflow** generate_stateflow_system6(int num, int maxState, double tUtil,int period_choice);
	/// if period_choice == 1, we use the function setup_periods_for_one_stateflow
	/// else if period_choice == 2, we use the function setup_periods_for_one_stateflow2
	/// else if period_choice == 3, we use the function setup_periods_for_one_stateflow3
	static Stateflow** generate_stateflow_system7(int num, int maxState, double tUtil, int period_choice);

	/// Corresponding to the generate_stateflow_system6, we add one parameter scc_probability
	static Stateflow** generate_stateflow_system8(int num, int maxState, double tUtil,int period_choice, double scc_probability);

	static Stateflow** generate_stateflow_system_for_exact_analysis(int num, int maxState, double tUtil,int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis(int num, int maxState, double tUtil, int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis2(int num, int maxState, double tUtil, int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis3(int num, int numState, double tUtil, int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis4(int num, int maxState, double tUtil, int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis5(int num, int maxState, double tUtil, int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis6(int num, int maxState, double tUtil, int period_choice, double scc_probability);
	static Stateflow** generate_stateflow_system_for_approximate_analysis7(int num, int maxState, double tUtil, int period_choice, double scc_probability);

	static Stateflow** generate_stateflow_system_for_approximate_analysis_with_periodicity(int num, int maxState, double tUtil, int period_choice, double scc_probability);
	

	/// Calculate the number of edges such taht the average degree for a vertex is around 3
	static int calculate_num_edge(int numNode);
};

#endif