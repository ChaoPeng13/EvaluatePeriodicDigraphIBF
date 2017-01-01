/* \file DigraphExperiments.h
 * This file contains the experiments for the digraph paper.
 * Two main parts: the improvment on ibf calculation and effect of digraph linear period
 * \author Chao Peng
 *
 * Changes
 *-------------------------------------------------------------------------------------------
 * 14-march-2016 : initial revision (CP)
 *
 */

#ifndef DIGRAPHEXPERIMENTS_H_
#define DIGRAPHEXPERIMENTS_H_

#include <direct.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

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

enum DigraphMethodChoice {
	DigraphExact,
	DigraphIBFGConstrainedDeadline, // graph searching, used in the IBF paper
	DigraphIBFLConstrainedDeadline, // linear property, used in the IBF paper
	DigraphIBFGLMADDeadline, // graph searching, used in the IBF paper
	DigraphIBFLLMADDeadline, // linear property, used in the IBF paper
	DigraphIBFGArbitraryDeadline, // graph searching, used in the IBF paper
	DigraphIBFLArbitraryDeadline, // linear property, used in the IBF paper
	DigraphRTIBFGConstrainedDeadline, // graph searching, used in the IBF paper
	DigraphRTIBFLConstrainedDeadline, // linear property, used in the IBF paper
	DigraphRBFG, // graph searching, used in the IBF paper
	DigraphRBFL, // linear property, used in the IBF paper
};

class DigraphExperiments {
public:
	static string digraphFuncNames[11];

	static bool doDigraphSchedAnalysis(DigraphMethodChoice dmc, Digraph** digraphs, int num);
	/// \brief the test system similar to the RTSS2014 by GUan et al.
	static void generateDigraphsForTestingAcceptanceRatio(string directory, int minUtil, int maxUtil, int stepUtil, int totalRun);
	static void testAcceptanceRatio(DigraphMethodChoice dmc, string directory, string file, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun);

	/// THis function generates the digraph systems to verify the improvment on ibf calculation
	/// But it's not used in our paper, because it needs a lot of time and the results are not smooth
	/// The digraph system will introduce too much randomness.
	static void generateDigraphsForTestingImprovementOnIBFCalculation(string directory);

	/// Similar to the above function, but here we use RandomGenerator::generate_mixed_digraph_system.
	static void generateDigraphsForTestingImprovementOnIBFCalculation(string directory,  int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int totalRun);
	/// This function implements the test for improvment on ibf calculation in digraph systems
	static void testImprovementOnIBFCalculationInDigraphSystem(string directory, string file, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int totalRun);

	/// We generate the random digraph system like the RTSJ2015 paper
	static void generateDigraphsForTestingImprovementOnIBFCalculation2(string directory, bool deadlineProperty, int maxNode,  int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun);
	static void generateDigraphsForTestingImprovementOnIBFCalculation2(string directory, bool deadlineProperty, int minNode, int maxNode,  int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun);
	static void generateDigraphsForTestingImprovementOnIBFCalculation2WithoutPointer(string directory, bool deadlineProperty, int minNode, int maxNode,  int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun);
	static void schedAnalysis(DigraphMethodChoice mc, string file, string directory,int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun);
	

	/// This function generate the digraph systems to verify the improvement of t_f
	/// But it's also not used sin our paper, because there is only 2% improvment (and straight)
	static void generateDigraphsForTestingImprovementOnTF(string directory, int minUtil, int maxUtil, int stepUtil, int totalRun);
	/// This function implements the test for improvment on t_f.
	/// Note that t_f only can be calculated in the system
	static void testImprovementOnTF(string directory, string file, int minUtil, int maxUtil, int stepUtil, int totalRun);

	/// This function generates only one digraph for the test of improvement on ibf calculation
	static void generateOneDigraphForTestingImprovemntOnIBFCalculation(string directory,int periodType, int maxRun);
	/// This function generates only one digraph for the test of improvement on ibf calculation, RTSJ16b
	static void generateOneDigraphForTestingImprovemntOnIBFCalculation2(string directory, int  utilInteger, int numNode, int maxRun);
	/// This function generates only one digraph for the test of improvement on ibf calculation, RTSJ16b
	static void generateOneDigraphForTestingImprovemntOnIBFCalculation2(string directory, int  utilInteger, int minNode, int maxNode, int maxRun);
	/// This function implements the test of improvement on ibf calculation in one digraph
	static void testImprovemntOnIBFCalculationInOneDigrah(string directory, string file, int periodType, int minRun, int maxRun, int stepRun);
	/// This function implements the test of improvement on ibf calculation in one digraph, RTSJ16b
	static void testImprovemntOnIBFCalculationInOneDigrah2(string directory, string file, int util, int numNode, int minRun, int maxRun, int stepRun);
	static void testLinearDefectForIBFRBF(string directory, string file, int periodType, int minRun, int maxRun, int stepRun);

	/// We implement this function to evaluate the effect of digraph task period, which is similar to Bini's paper
	/// There are two tasks, one is a digraph (with high priority) and another is a periodic task (with low priority).
	/// The digraph task is generated like Guan's paper
	/// T_h/T_l = [0,1]
	/// U_h/U_l = 0.25 and U=U_h+U_l={0.2,0.4,0.6,0.8}
	static void generateDigraphsForTestingEffectOfDigraphPeriod(string directory, int maxRun);
	/// We implement this function to evaluate the effect of digraph task period, which is similar to Bini's paper
	/// There are two tasks, one is a digraph (with high priority) and another is a periodic task (with low priority).
	/// The digraph task is generated like the RTSJ2015
	/// T_h/T_l = [0,1]
	/// U_h/U_l = 0.25 and U=U_h+U_l={0.2,0.4,0.6,0.8}
	static void generateDigraphsForTestingEffectOfDigraphPeriod2(string directory, int fixedNodeNum, int minRun, int maxRun);
	static void generateDigraphsForTestingEffectOfDigraphPeriodUtil(string directory, int fixedNodeNum, int minUtil, int maxUtil, int stepUtil, int minRun, int maxRun);
	static void testEffectOfDigraphPeriod(string directory, string file, int startRun, int endRun, int choice);
	static void testEffectOfDigraphPeriodUtil(string directory, string file, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun);
	static bool testEffectOfDigraphPeriod(vector<double> &IBFResult, vector<double> &RBFResults, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double>& ratios,  const char *p, int util, double &exactUtil); 
	static void testEffectOfDigraphPeriod2(vector<double> &ExactResult, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double>& ratios,  const char *p, int util, double &exactUtil);
	static bool testEffectOfDigraphPeriodWithHarmonicSet(vector<double> &IBFResult, vector<double> &RBFResults, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double>& ratios,  const char *p, int util, double &exactUtil); 
	static bool testEffectOfDigraphPeriodWithHarmonicSetUtil(vector<double> &IBFResult, vector<double> &RBFResults, vector<double> &LUBRBFResult, vector<double> &LUBIBFResult, set<double> ratios,  const char *p, int util, double &exactUtil);
	
	/// We implement this function to evaluate the effect of task number, which is similar to Bini's paper
	/// There are n tasks, n-1 digraphs (whih high priority) and the remaining one is a periodic task (with low priority)
	/// T_h = \max{T_1,\ldots,T_{n-1}}
	static void generateDigraphsForTestingEffectOfTaskParameters(string directory, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int maxRun);
	static void testEffectOfTaskParameters(string directory, string file, bool isExact, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun, int choice);
	/// We implement this function to evaluate the effect of task number, which is similar to Bini's paper
	/// There are n+1 tasks, n digraphs (whih high priority) and the remaining one is a periodic task (with low priority)
	/// T_h = \max{T_1,\ldots,T_{n}}
	static void generateDigraphsForTestingEffectOfTaskParameters2(string directory,int maxNodeNum, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun);
	static void testEffectOfTaskParameters2(string directory, string file, double factor, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun, int choice);
	static void testEffectOfTaskParameters2Util(string directory, string file, double factor, int minNum, int maxNum, int stepNum, int minUtil, int maxUtil, int stepUtil, int startRun, int endRun, int choice);


	static void outputResults(ofstream& fout, vector<vector<double>> results, string name);
	static void outputResults2(ofstream& fout, vector<vector<double>> results, string name);
	static void outputResults(ofstream& fout, map<int,map<int,double>> results, string name);
	static void outputResults(ofstream& fout, map<int,map<int,int>> results, string name);
	static void outputResults(ofstream& fout, map<int,double> results, string name);
	static void outputResults(ofstream& fout, map<int,int> results, string name);
};

#endif