/* \file FSMExperiments.h
 * This file contains the experiments for the FSM papers.
 * The main parts are: the schedubility analyses including 
 * ->RF/IF-SO 
 * ->RF/IF-AO
 * ->RBF-SO 
 * ->IBF-SO 
 * ->RBF-AO1 
 * ->RBF-AO2 
 * ->LUBRBF 
 * ->LUBIBF 
 * \author Chao Peng
 *
 * Changes
 *-------------------------------------------------------------------------------------------
 * 14-march-2016 : initial revision (CP)
 *
 */

#ifndef FSMEXPERIMENTS_H_
#define FSMEXPERIMENTS_H_

#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include <future>
#include <chrono>

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

typedef bool (*PointFunc) (Stateflow** sfs, int n, int choice);

enum MethodChoice {
	EXACTSO = 0, // Exact schedulability analysis on FSMs with static offsets
	EXACTAO1, // Exact schedulability analysis on (precise) digraphs transformed with action instances
	EXACTAO2, // Exact schedulability analysis on (simple) digraphs transfromed with actions
	IBFSO, // IBF schedulability analysis on FSMs with static offsets
	RBFSO, // RBF schedulability analysis on FSMs with static offsets
	IBFSOP, // IBF schedulability analysis on FSMs with static offsets using the linear periodicity
	RBFSOP, // RBF schedulability analysis on FSMs with static offsets using the linear periodicity
	IBFAOFSM, // IBF schedulability analysis with arbitrary offsets on FSMs
	RBFAOFSM, // RBF schedulability analysis with arbitrary offsets on FSMs
	IBFAOFSMDP, // IBF schedulability analysis with arbitrary offsets on FSMs using dynamic programming techniques
	RBFAOFSMDP, // RBF schedulability analysis with arbitrary offsets on FSMs using dynamic programming techniques
	IBFAO1G, // IBF schedulability analysis on (precise) digraphs using Guan's algorithm
	RBFAO1G, // RBF schedulability analysis on (precise) digraphs using Guan's algorithm
	IBFAO1DP, // IBF schedulability analysis on (precise) digraphs using dynamic programming techniques
	RBFAO1DP, // RBF schedulability analysis on (precise) digraphs using dynamic programming techniques
	IBFAO1L, // IBF schedulability analysis on (precise) digraphs using the linearization
	RBFAO1L, // RBF schedulability analysis on (precise) digraphs using the linearization
	IBFAO1LSCC, // IBF schedulability analysis on (strongly-connected precise) digraphs using the linearization
	RBFAO1LSCC, // RBF schedulability analysis on (strongly-connected precise) digraphs using the linearization
	IBFAO2G, // IBF schedulability analysis on (simple) digraphs using Guan's algorithm
	RBFAO2G, // RBF schedulability analysis on (simple) digraphs using Guan's algorithm
	IBFAO2DP, // IBF schedulability analysis on (simple) digraphs using the dynamic programming techniques
	RBFAO2DP, // RBF schedulability analysis on (simple) digraphs using the dynamic programming techniques
	IBFAO2L, // IBF schedulability analysis on (simple) digraphs using the linearization
	RBFAO2L, // RBF schedulability analysis on (simple) digraphs using the linearization
	LUBIBF, // LUBIBF schedulability analysis on FSMs with static offsets
	LUBRBF, // LUBRBF schedulability analysis on FSMs with static offsets
	LUBCSUM, // LUBCSUM schedulability analysis on FSMs with static offsets
};

class FSMExperiments {
public:
	static string funcNames[28];
	static void generateRandomSystemsForExactAnalysis(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int minRun, int maxRun);
	static void generateRandomSystemsForApproximateAnalysis(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int maxRun);
	static void generateRandomSystemsForApproximateAnalysisWithPeriodicity(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int maxRun);
	static void generateRandomSystemsForApproximateAnalysisWithPeriodicity2(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int minRun, int maxRun);
	/// test all the analysis methods
	/// We use the parameter myProperty to separate two types of analyses
	/// myProperty = true, do exact analysis
	/// otherwise, do approximate analysis
	static void SchedAnalysis(MethodChoice mc, string file, string directory,int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int startRun, int endRun);
	
	static void generateRandomSystemsForWeightedSchedAnalysis(string directory, int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int minScc, int maxScc, int stepScc, int maxRun); 
	static void WeightedSchedulabilityAnalysis(MethodChoice mc,string directory, string file, int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int minScc, int maxScc, int stepScc, int maxRun);

	static void generateRandomSystemsForScalabilityTest(string directory,int period_choice, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int minScc, int maxScc, int stepScc, int minState, int maxState, int stepState, int minRun, int maxRun);
	/// Every time, We try to test the scalability for only one analysis method 
	static void SchedAnalysisForScalabilityTest(MethodChoice mc, bool myProperty, string file, string directory, int choice, int num, int util, int scc, int minState, int maxState, int stepState, int startRun, int endRun);

	static bool doSchedAnalysis(MethodChoice mc, Stateflow** sfs, int num, ostream& out);
	static void output(map<int,map<int,vector<double>>> allRecTimes, string name, ostream& out);
};

#endif