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

int main(int argc, char* argv[])
{
#ifdef GOOGLETEST
	// Run Google test
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
#endif
	string index = "20161220";
	string result = "Results\\";
	string resultDir = "Results";
	string suffix = ".res";
	int times = 2;
	bool deadlineProperty = true;

	// mkdir result directory
	if ( _mkdir(resultDir.c_str()) == 0) {
		cout<<"Directory: "<<resultDir<<" was successfully created."<<endl;
	}

	//=====================================================================================
	// Calculation Efficiency Improvement from Linear Periodicity
	//=====================================================================================

#if 1 // fixed vertices number, {5,10,15,20,25}
	for (int numNode = 5; numNode <= 25; numNode += 5) {
		DigraphExperiments::generateOneDigraphForTestingImprovemntOnIBFCalculation2("Digraphs"+index+"Node"+Utility::int_to_string(numNode),50,numNode,times);
		DigraphExperiments::testImprovemntOnIBFCalculationInOneDigrah2("Digraphs"+index+"Node"+Utility::int_to_string(numNode),result+"DigraphTimeCalIBF"+index+"Node"+Utility::int_to_string(numNode),50,numNode,0,times-1,1);
	}
#endif

#if 1 // random vertices number between 10 and 25
	DigraphExperiments::generateOneDigraphForTestingImprovemntOnIBFCalculation2("Digraphs"+index+"NodeRandom",50,10,25,times);
	DigraphExperiments::testImprovemntOnIBFCalculationInOneDigrah2("Digraphs"+index+"NodeRandom",result+"DigraphTimeCalIBF"+index+"NodeRandom",50,1025,0,times-1,1);
#endif

	//=====================================================================================
	// Analysis Efficiency Improvement from Linear Periodicity
	//=====================================================================================

#if 1 // Evaluate the analysis efficiency versus the total utilization
	DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2("ImprovementUtil"+index,deadlineProperty, 10, 25,20,20,1,50,95,5,0,times-1);
	DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"ImprovementUtil"+index,"ImprovementUtil"+index,20,20,1,50,95,5,0,times-1);
	SchedulabilityAnalysis::clear_all();
	DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"ImprovementUtil"+index,"ImprovementUtil"+index,20,20,1,50,95,5,0,times-1);
	SchedulabilityAnalysis::clear_all();
#endif


#if 1 // Evaluate the analysis efficiency versus the task number
	for (int fixedUtil = 90; fixedUtil <=90; fixedUtil += 10) {
		DigraphExperiments::generateDigraphsForTestingImprovementOnIBFCalculation2("ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),deadlineProperty,10,25,5,40,5,fixedUtil,fixedUtil,5,0,times-1);
		DigraphExperiments::schedAnalysis(DigraphIBFGArbitraryDeadline, result+"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),5,40,5,fixedUtil,fixedUtil,5,0,times-1);
		SchedulabilityAnalysis::clear_all();
		DigraphExperiments::schedAnalysis(DigraphIBFLArbitraryDeadline, result+"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),"ImprovementNum"+index+"Util"+Utility::int_to_string(fixedUtil),5,40,5,fixedUtil,fixedUtil,5,0,times-1);
		SchedulabilityAnalysis::clear_all();
	}
#endif

	//=====================================================================================
	// Effect of Task Periods with Two-tasks Systems
	//=====================================================================================

#if 1 // versus the ratio $T_h\T_l$ given the fixing utilization and vertex number
	int fixedNodeNumbers[5] = {2,5,10,20,50};
	for (int i = 0; i < 5; i ++) {
		int fixedNodeNum = fixedNodeNumbers[i];
		DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriod2("DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum),fixedNodeNum,0,times-1);
		DigraphExperiments::testEffectOfDigraphPeriod("DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum), result+"DigraphPeriod"+index+"Node"+Utility::int_to_string(fixedNodeNum),0,times-1,3);
	}
#endif

#if 1 // versus the utilization given the fixing $T_h\T_l$ and vertex number
	int fixedNodeNumbers2[2] = {2,20};
	times = 1000;
	for (int i = 0; i < 2; i ++) {
		int fixedNodeNum = fixedNodeNumbers2[i];
		DigraphExperiments::generateDigraphsForTestingEffectOfDigraphPeriodUtil("DigraphPeriodUtil"+index+"Node"+Utility::int_to_string(fixedNodeNum),fixedNodeNum,10,90,10,0,times-1);
		DigraphExperiments::testEffectOfDigraphPeriodUtil("DigraphPeriodUtil"+index+"Node"+Utility::int_to_string(fixedNodeNum), result+"DigraphPeriodUtil"+index+"Node"+Utility::int_to_string(fixedNodeNum),10,90,10,0,times-1);
	}
#endif

	//=====================================================================================
	// Effect of Task System Parameters with Larger Systems
	//=====================================================================================

#if 1 // versus the task number given the fixing utilization and factor
	int utils[3] = {30,50,70}; 
	int factors[3] = {20,50,80};

	for (int i=0; i<3; i++) {
		int util = utils[i];
		DigraphExperiments::generateDigraphsForTestingEffectOfTaskParameters2("DigraphParameters"+index,40,5,40,5,util,util,10,0,times-1);
		for (int j=0; j<3; j++) {
			int factor = factors[j];
			DigraphExperiments::testEffectOfTaskParameters2("DigraphParameters"+index, result+"DigraphParameters"+index+"Util"+Utility::int_to_string(util)+"factor"+Utility::int_to_string(factor),0.01*factor,5,40,5,util,util,10,0,times-1,1);
		}
	}
#endif

	return 0;
}