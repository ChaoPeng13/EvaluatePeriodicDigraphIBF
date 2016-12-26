/* \file SchedulabilityAnalysis.h
 * this file implements the schedulability analysis of digraphs and statflows
 * \author Chao Peng
 *
 * Changes
 * -----------
 * 07-sept-2015 : initial revision (CP)
 *
 */

#ifndef SCHEDULABILITYANALYSIS_H_
#define SCHEDULABILITYANALYSIS_H_

//#include <vld.h>

#include "Digraph.h"
#include "Stateflow.h"

class SchedulabilityAnalysis {
public:
	//===================================================================================
	// Digraph information
	//===================================================================================
	// Statistics Data
	static double tfRatio0;
	static double tfRatio1;

	static int nDigraphExact;
	static int nDigraphIBF;
	static int nDigraphRBF;

	static int nIrreducibleDigraphs;
	static int nDigraphNode;
	static int nDigraphEdge;
	// Statistics Time
	static double tDigraphCalLinearFactor;
	static double tDigraphCalLinearPeriod;
	static double tDigraphCalLinearDefect;
	static double tDigraphCalSum;
	static double tDigraphCalLinearBounds;
	static double tDigraphCalTF0;
	static double tDigraphCalTF1;
	static double tDigraphCalTF2;

	static double tDigraphExact;
	static double tDigraphRBF;
	static double tDigraphIBF;
	// Vectors
	static vector<double> vec_tfRatio0;
	static vector<double> vec_tfRatio1;

	static vector<int> vec_nDigraphExact;
	static vector<int> vec_nDigraphIBF;
	static vector<int> vec_nDigraphRBF;

	static vector<int> vec_nIrreducibleDigraphs;
	static vector<int> vec_nDigraphNode;
	static vector<int> vec_nDigraphEdge;

	static vector<double> vec_tDigraphCalLinearFactor;
	static vector<double> vec_tDigraphCalLinearPeriod;
	static vector<double> vec_tDigraphCalLinearDefect;
	static vector<double> vec_tDigraphCalSum;
	static vector<double> vec_tDigraphCalLinearBounds;
	static vector<double> vec_tDigraphCalTF0;
	static vector<double> vec_tDigraphCalTF1;
	static vector<double> vec_tDigraphCalTF2;

	static vector<double> vec_tDigraphExact;
	static vector<double> vec_tDigraphRBF;
	static vector<double> vec_tDigraphIBF;
	//===================================================================================
	// Stateflow information
	//===================================================================================
	static int nStateflows;
	static double totalUtilization;
	static double totalUtilization2;
	static double avgDegree;
	// statistics Results
	// with static offsets
	static int nTimeout;
	static int nExactStaticOffset;
	static int nRBFStaticOffset;
	static int nIBFStaticOffset;
	static int nRBFArbitraryOffset;
	static int nIBFArbitraryOffset;

	// with arbitrary offsets
	static int nExactArbitraryOffsetBySimpleDigraph; // Exact schedulability analysis on (simple) digraphs transformed with action instances
	static int nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static int nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static int nRBFArbitraryOffsetBySimpleDigraphWithLinearization;
	static int nIBFArbitraryOffsetBySimpleDigraphWithLinearization;

	static int nExactArbitraryOffsetByPreciseDigraph; // Exact schedulability analysis on (precise) digraphs transformed with action instances
	static int nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static int nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static int nRBFArbitraryOffsetByPreciseDigraphWithLinearization;
	static int nIBFArbitraryOffsetByPreciseDigraphWithLinearization;

	static int nLinearUpperRBF;
	static int nLinearUpperIBF;
	static int nLinearUpperCSUM;
	// Statistics Data
	static double bpRatio0; // c^rbf+c^dbf /c^sum 
	static double bpRatio1; // c^ibf+c^dbf / c^sum

	static int nIrreducibleStateflows;
	// Statistics Time
	static double tPrepareStateflow; // time for preparing stateflow such as calculating gcd and hyperperiod, generating exeuction request matrix and etc. 

	static double tCalCSum; // time for calculating c^sum
	static double tCalLinearBounds; // time for calculating the linear upper bounds
	static double tCalTF0; // time for calculating t_f baed on c^sum
	static double tCalTF1; // time for calculating t_f based on c^rbf and c^dbf
	static double tCalTF2; // time for calculating t_f based on c^ibf and c^dbf
	static double tCalNonPeriod; // time for calculating the matrix power without periodicity property
	static double tCalPeriod; // time for calculating the matrix power with periodicity property
	static double tCalDiffPeriod; // the different time for calculating the matrix power
	static double tCalRequestExecutionMatrixWithoutPeriodicityProperty;
	static double tCalRequestExecutionMatrixWithPeriodicityProperty;
	static double tGenerateCriticalActionPairs;
	static double tGenerateRequestFunctions;
	static double tGenerateRequestFunctionAbstractTree;
	
	// with static offsets
	static double tExactStaticOffset;
	static double tRBFStaticOffset;
	static double tIBFStaticOffset;
	static double tRBFArbitraryOffset;
	static double tIBFArbitraryOffset;

	// with arbitrary offsets
	static double tExactArbitraryOffsetBySimpleDigraph; // Exact schedulability analysis on (simple) digraphs transformed with action instances
	static double tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static double tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static double tRBFArbitraryOffsetBySimpleDigraphWithLinearization;
	static double tIBFArbitraryOffsetBySimpleDigraphWithLinearization;

	static double tExactArbitraryOffsetByPreciseDigraph; // Exact schedulability analysis on (precise) digraphs transformed with action instances
	static double tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static double tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static double tRBFArbitraryOffsetByPreciseDigraphWithLinearization;
	static double tIBFArbitraryOffsetByPreciseDigraphWithLinearization;

	static double tLinearUpperRBF;
	static double tLinearUpperIBF;
	static double tLinearUpperCSUM;
	
	// Vectors
	static vector<int> vec_nStateflows;
	static vector<double> vec_totalUtilization;
	static vector<double> vec_totalUtilization2;
	static vector<double> vec_avgDegree;

	static vector<int> vec_nTimeout;
	static vector<int> vec_nExactStaticOffset;
	static vector<int> vec_nRBFStaticOffset;
	static vector<int> vec_nIBFStaticOffset;
	static vector<int> vec_nRBFArbitraryOffset;
	static vector<int> vec_nIBFArbitraryOffset;
	static vector<int> vec_nExactArbitraryOffsetBySimpleDigraph; // Exact schedulability analysis on (simple) digraphs transformed with action instances
	static vector<int> vec_nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static vector<int> vec_nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static vector<int> vec_nRBFArbitraryOffsetBySimpleDigraphWithLinearization;
	static vector<int> vec_nIBFArbitraryOffsetBySimpleDigraphWithLinearization;
	static vector<int> vec_nExactArbitraryOffsetByPreciseDigraph; // Exact schedulability analysis on (precise) digraphs transformed with action instances
	static vector<int> vec_nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static vector<int> vec_nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static vector<int> vec_nRBFArbitraryOffsetByPreciseDigraphWithLinearization;
	static vector<int> vec_nIBFArbitraryOffsetByPreciseDigraphWithLinearization;
	static vector<int> vec_nLinearUpperRBF;
	static vector<int> vec_nLinearUpperIBF;
	static vector<int> vec_nLinearUpperCSUM;

	static vector<double> vec_bpRatio0;
	static vector<double> vec_bpRatio1;
	static vector<int> vec_nIrreducibleStateflows;

	static vector<double> vec_tPrepareStateflow;

	static vector<double> vec_tCalCSum;
	static vector<double> vec_tCalLinearBounds;
	static vector<double> vec_tCalTF0;
	static vector<double> vec_tCalTF1;
	static vector<double> vec_tCalTF2;
	static vector<double> vec_tCalNonPeriod;
	static vector<double> vec_tCalPeriod;
	static vector<double> vec_tCalDiffPeriod;
	static vector<double> vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty;
	static vector<double> vec_tCalRequestExecutionMatrixWithPeriodicityProperty;
	static vector<double> vec_tGenerateCriticalActionPairs;
	static vector<double> vec_tGenerateRequestFunctions;
	static vector<double> vec_tGenerateRequestFunctionAbstractTree;
	
	static vector<double> vec_tExactStaticOffset;
	static vector<double> vec_tRBFStaticOffset;
	static vector<double> vec_tIBFStaticOffset;
	static vector<double> vec_tRBFArbitraryOffset;
	static vector<double> vec_tIBFArbitraryOffset;
	static vector<double> vec_tExactArbitraryOffsetBySimpleDigraph; // Exact schedulability analysis on (simple) digraphs transformed with action instances
	static vector<double> vec_tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static vector<double> vec_tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
	static vector<double> vec_tRBFArbitraryOffsetBySimpleDigraphWithLinearization;
	static vector<double> vec_tIBFArbitraryOffsetBySimpleDigraphWithLinearization;
	static vector<double> vec_tExactArbitraryOffsetByPreciseDigraph; // Exact schedulability analysis on (precise) digraphs transformed with action instances
	static vector<double> vec_tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static vector<double> vec_tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
	static vector<double> vec_tRBFArbitraryOffsetByPreciseDigraphWithLinearization;
	static vector<double> vec_tIBFArbitraryOffsetByPreciseDigraphWithLinearization;
	static vector<double> vec_tLinearUpperRBF;
	static vector<double> vec_tLinearUpperIBF;
	static vector<double> vec_tLinearUpperCSUM;
                      
	//=====================================================================================
	// Method start
	//=====================================================================================
	static void save_results();
	static void set_zero();
	static void clear_all();
	static void reset();
	static void output_one_vector(ofstream& fout, string sVec, vector<double> vec);
	static void output_one_vector(ofstream& fout, string sVec, vector<int> vec);
	static void output_vectors(ofstream& fout);
	/* Choice used to identify which busy length to choose 
	 * 0 -> original busy length, i.e. csum
	 * 1 -> C^rbf and C^dbf
	 * 2 -> C^ibf and C^dbf
	 */
	
	// prepare t_f 
	static void prepare_all_digraphs(Digraph** digraphs, int n, int choice);
	static void prepare_all_digraphs2(Digraph** digraphs, int n); 

	static bool isAbstract(vector<DigraphRFNode*> arfs);

	static void generate_critical_vertices(Digraph** digraphs, int n);
	static void generate_critical_vertices_for_reachability_digraphs(Digraph** digraphs, int n);
	
	static double calculate_digraph_exact_response_time(Digraph** digraphs, int i, int wcet);
	static double calculate_digraph_exact_response_time(vector<DigraphRFNode*> combination, int wcet);

	static bool digraph_exact_analysis_timeout(Digraph** digraphs, int n, int choice);
	static bool digraph_exact_analysis(Digraph** digraphs, int n, int choice);
	static bool digraph_exact_analysis(Digraph** digraphs, int i);
	static bool digraph_exact_analysis(vector<DigraphRFNode*> combination, int wcet, int deadline);
	//static bool digraph_exact_analysis(vector<DigraphRFNode*> roots, int wcet, int deadline);
	//static bool digraph_exact_analysis(vector<vector<DigraphRFNode*>> combinations, int wcet, int deadline);
	

	/* linearProperty used to identify how to calculate rbf and ibf
	 * true -> with execution matrix
	 * false -> without execution matrix
	 */
	static bool digraph_rbf_analysis(Digraph** digraphs, int n, int choice, bool linearization);
	static bool digraph_rbf_analysis(Digraph** digraphs, int i, bool linearization);
	static bool digraph_rbf_analysis(Digraph** digraphs, int i, bool linearization, int wcet, int deadline);

	static bool digraph_rbf_analysis2(Digraph** digraphs, int n, int choice);
	static bool digraph_rbf_analysis2(Digraph** digraphs, int i);
	static bool digraph_rbf_analysis2(Digraph** digraphs, int i, int wcet, int deadline);

	static bool digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice,bool periodicityProperty,bool deadlineProperty);
	static bool digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i);
	static bool digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline);

	static bool reachability_digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice);
	static bool reachability_digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i);
	static bool reachability_digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline);

	static bool digraph_ibf_analysis(Digraph** digraphs, int n, int choice, bool linearization);
	static bool digraph_ibf_analysis(Digraph** digraphs, int i, bool linearization);
	static bool digraph_ibf_analysis(Digraph** digraphs, int i, bool linearization, int wcet, int deadline);

	static bool digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice,bool periodicityProperty,bool deadlineProperty);
	static bool digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i, bool deadlineProperty);
	static bool digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline);

	static bool digraph_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int n, int choice, bool periodicityProperty);
	static bool digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i);
	static bool digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i, double maxL, Node* node);
	static double digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i, double maxL, RequestTriple* rt);
	static bool digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i, double maxL, int wcet, int deadline, int relaseTime);

	static bool digraph_rt_ibf_analysis_dynamic_programming_under_constrained_deadline(Digraph** digraphs, int n, int choice, bool periodicityProperty);
	static double digraph_rt_ibf_analysis_dynamic_programming_under_constrained_deadline(Digraph** digraphs, int i, int wcet, double maxRT);

	static bool reachability_digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice);
	static bool reachability_digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i);
	static bool reachability_digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline);

	static bool digraph_ibf_analysis_with_rbf_linear_defect(Digraph** digraphs, int n, int choice);
	static bool digraph_ibf_analysis_with_rbf_linear_defect(Digraph** digraphs, int i);
	static bool digraph_ibf_analysis_with_rbf_linear_defect(Digraph** digraphs, int i, int wcet, int deadline);

	/* Choice used to identify which busy length to choose 
	 * 0 -> original busy length, i.e. csum
	 * 1 -> C^rbf and C^dbf
	 * 2 -> C^ibf and C^dbf
	 */

	/// \brief prepare linear upper bounds and the maximal time length
	static void prepare_all_stateflows(Stateflow** sfs, int n, int choice);
	/// \brief prepare rbf or ibf and tf is calculated by the deadlines not linear upper bounds
	static void prepare_all_stateflows2(Stateflow** sfs, int n, bool RBF, bool periodicity);

	// reset some containters for all the stateflows
	static void reset_calculating_containers(Stateflow** sfs, int n);

	// calculate the request execution matrices, in particular to statistic the runtimes
	static void calculate_request_execution_matrices_without_periodicity_property(Stateflow** sfs, int n, int choice);
	static void calculate_request_execution_matrices_with_periodicity_property(Stateflow** sfs, int n, int choice);

	static void generate_critical_action_pair(Stateflow** sfs, int n);

	static bool generate_request_function(Stateflow** sfs, int n, bool output);
	static bool generate_request_function(Stateflow* sf, int i, int maxDeadline, bool output);

	static void generate_request_function_for_each_action_instance(Stateflow** sfs, int n, bool output);
	static void generate_request_function_for_each_action_instance(Stateflow** sfs, int i, int tk, int maxDeadline, bool output);
	static void generate_request_function_for_each_action_instance(Stateflow* sf, int i, int tk, int maxDeadline, bool output);

	static void intra_loop_dominance(vector<RequestFunction>& vec_rf);
	static void end_loop_dominance(vector<RequestFunction>& vec_rf);

	static void generate_request_function_abstract_tree(Stateflow** sfs, int n, bool output);
	static void generate_request_function_abstract_tree(Stateflow** sf, int i, int stime, int deadline, bool output);
	static void generate_request_function_abstract_tree(Stateflow* sf, int stime, int deadline, bool output);

	static void generate_request_function_abstract_tree_for_each_action_instance(Stateflow** sfs, int n, bool output);
	static void generate_request_function_abstract_tree_for_each_action_instance(Stateflow* sf);


	static void generate_one_request_function_abstract_tree(Stateflow** sfs, int n, bool output);
	static void generate_one_request_function_abstract_tree(Stateflow* sf);

	static bool isAbstract(vector<RFNode*> vec_rfnode);

	// Exact schedulability analysis
	static bool exact_sched_analysis(Stateflow** sfs, int n, int choice,bool output, ostream& out);
	static bool exact_sched_analysis(Stateflow** sfs, int i, int s, vector<ActionPair> vec_ap);
	static bool exact_sched_analysis(Stateflow** sfs, int i, int s, ActionPair ap);
	static bool exact_sched_analysis(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int stime, int deadline, int gcd);
	static bool exact_sched_analysis(vector<vector<RFNode*>> Omega, int wcet, int stime, int deadline, int gcd);
	static bool exact_sched_analysis(vector<RFNode*> vec_rfnode, int wcet, int stime, int deadline, int gcd);
	static int calculate_remaining_workload(vector<RFNode*> vec_rfnode, int stime, int gcd);

	// Exact analysis with arbitrary offset
	static bool exact_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static bool exact_analysis_arbitrary_offset(Stateflow** sfs, int i, vector<ActionPair> vec_ap);
	static bool exact_analysis_arbitrary_offset(Stateflow** sfs, int i, ActionPair ap);
	static bool exact_analysis_arbitrary_offset(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int deadline, int gcd);
	static bool exact_analysis_arbitrary_offset(vector<vector<RFNode*>> Omega, int wcet, int deadline, int gcd);
	static bool exact_analysis_arbitrary_offset(vector<RFNode*> vec_rfnode, int wcet, int deadline, int gcd);

	// Exact analysis on simple digraph models
	static bool exact_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice);
	// Exact analysis on precise digraph models
	static bool exact_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice);
	static bool exact_analysis_arbitrary_offset_by_precise_digraphs_timeout(Stateflow** sfs, int n, int choice, int timeout);

	// Approximate schedulability analysis based on RBF with static offsets
	static bool rbf_analysis_static_offset(Stateflow** sfs, int n, int choice, bool periodicity);
	static bool rbf_analysis_static_offset_index(Stateflow** sfs, int i);
	static bool rbf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);
	static double calculate_remaining_rbf_workload(Stateflow** sfs, int i, int stime, int gcd);

	// Approximate schedulability analysis based on IBF with static offsets
	static bool ibf_analysis_static_offset(Stateflow** sfs, int n, int choice, bool periodicity);
	static bool ibf_analysis_static_offset_index(Stateflow** sfs, int i);
	static bool ibf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);
	static double calculate_remaining_ibf_workload(Stateflow** sfs, int i, int stime, int gcd);

	// Approximate schedulability analysis based on RBF with arbitrary offsets
	static bool rbf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static bool rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i);
	static bool rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline);

	// RBF analysis on simple models
	static bool rbf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice,bool linearization);
	static bool rbf_analysis_arbitrary_offset_by_simple_digraphs_DP(Stateflow** sfs, int n, int choice, bool periodicityProperty);
	// RBF analysis on precise models
	static bool rbf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice, bool linearization);
	static bool rbf_analysis_arbitrary_offset_by_precise_digraphs_DP(Stateflow** sfs, int n, int choice, bool periodicityProperty);
	static bool rbf_analysis_arbitrary_offset_by_strongly_connected_precise_digraphs_DP(Stateflow** sfs, int n, int choice, bool periodicityProperty);
	static bool rbf_analysis_arbitrary_offset_by_reachability_digraphs_DP(Stateflow** sfs, int n, int choice, bool periodicityProperty);
	
	// Approximate schedulability analysis based on IBF with arbitrary offsets
	static bool ibf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice);
	static bool ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i);
	static bool ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline);

	// IBF analysis on simple models
	static bool ibf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice, bool linearization);
	static bool ibf_analysis_arbitrary_offset_by_simple_digraphs_DP(Stateflow** sfs, int n, int choice, bool periodicityProperty);
	// IBF analysis on precise models
	static bool ibf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice, bool linearization);
	static bool ibf_analysis_arbitrary_offset_by_precise_digraphs_DP(Stateflow** sfs, int n, int choicebool, bool periodicityProperty);
	static bool ibf_analysis_arbitrary_offset_by_strongly_connected_precise_digraphs_DP(Stateflow** sfs, int n, int choicebool, bool periodicityProperty);
	static bool ibf_analysis_arbitrary_offset_by_reachability_digraphs_DP(Stateflow** sfs, int n, int choice, bool periodicityProperty);

	// Linear schedulability analysis based on RBF with static offset
	static bool lu_rbf_sched_analysis(Stateflow** sfs, int n, int choice);
	static bool lu_rbf_sched_analysis_index(Stateflow** sfs, int i);
	static bool lu_rbf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);

	// Linear schedulability analysis based on IBF with static offset
	static bool lu_ibf_sched_analysis(Stateflow** sfs, int n, int choice);
	static bool lu_ibf_sched_analysis_index(Stateflow** sfs, int n);
	static bool lu_ibf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);

	// Linear schedulability analysis based on CSUM with static offset
	static bool lu_csum_sched_analysis(Stateflow** sfs, int n, int choice);
	static bool lu_csum_sched_analysis_index(Stateflow** sfs, int n);
	static bool lu_csum_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline);

	// output
	static void output_critical_action_pair(Stateflow** sfs, int n, ostream& out);

	static void output_critical_request_function(Stateflow** sfs, int n, ostream& out);
	static void output_critical_request_function(Stateflow* sf, ostream& out);
	static void output_critical_request_function(vector<RequestFunction> vec_rf, ostream& out);

	static void output_request_function_abstract_tree(Stateflow** sfs, int n, ostream& out);
	static void output_request_function_abstract_tree(Stateflow* sf, ostream& out);
	static void output_request_function_abstract_tree(AbstractRequestFunctionTree arft, ostream& out);
	
};

#endif