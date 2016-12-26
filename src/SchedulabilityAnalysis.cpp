#include "SchedulabilityAnalysis.h"
#include "Timer.h"

#include <algorithm>
#include <Windows.h>
#include <future>
#include <chrono>

double SchedulabilityAnalysis::tfRatio0 = 0;
double SchedulabilityAnalysis::tfRatio1 = 0;

int SchedulabilityAnalysis::nDigraphExact = 0;
int SchedulabilityAnalysis::nDigraphIBF = 0;
int SchedulabilityAnalysis::nDigraphRBF = 0;

int SchedulabilityAnalysis::nIrreducibleDigraphs =0;
int SchedulabilityAnalysis::nDigraphNode = 0;
int SchedulabilityAnalysis::nDigraphEdge = 0;

double SchedulabilityAnalysis::tDigraphCalLinearFactor = 0;
double SchedulabilityAnalysis::tDigraphCalLinearPeriod = 0;
double SchedulabilityAnalysis::tDigraphCalLinearDefect = 0;
double SchedulabilityAnalysis::tDigraphCalSum = 0;
double SchedulabilityAnalysis::tDigraphCalLinearBounds = 0;
double SchedulabilityAnalysis::tDigraphCalTF0 = 0;
double SchedulabilityAnalysis::tDigraphCalTF1 = 0;
double SchedulabilityAnalysis::tDigraphCalTF2 = 0;

double SchedulabilityAnalysis::tDigraphExact = 0;
double SchedulabilityAnalysis::tDigraphRBF = 0;
double SchedulabilityAnalysis::tDigraphIBF = 0;

int SchedulabilityAnalysis::nStateflows = 0;
double SchedulabilityAnalysis::totalUtilization = 0;
double SchedulabilityAnalysis::totalUtilization2 = 0;
double SchedulabilityAnalysis::avgDegree = 0;

int SchedulabilityAnalysis::nTimeout = 0;
int SchedulabilityAnalysis::nExactStaticOffset = 0;
int SchedulabilityAnalysis::nRBFStaticOffset = 0;
int SchedulabilityAnalysis::nIBFStaticOffset = 0;
int SchedulabilityAnalysis::nRBFArbitraryOffset = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffset = 0;
int SchedulabilityAnalysis::nExactArbitraryOffsetBySimpleDigraph = 0; // Exact schedulability analysis on (simple) digraphs transformed with action instances
int SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
int SchedulabilityAnalysis::nRBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
int SchedulabilityAnalysis::nExactArbitraryOffsetByPreciseDigraph = 0; // Exact schedulability analysis on (precise) digraphs transformed with action instances
int SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
int SchedulabilityAnalysis::nRBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
int SchedulabilityAnalysis::nIBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
int SchedulabilityAnalysis::nLinearUpperRBF = 0;
int SchedulabilityAnalysis::nLinearUpperIBF = 0;
int SchedulabilityAnalysis::nLinearUpperCSUM = 0;

double SchedulabilityAnalysis::bpRatio0 = 0;
double SchedulabilityAnalysis::bpRatio1 = 0;

int SchedulabilityAnalysis::nIrreducibleStateflows = 0;

double SchedulabilityAnalysis::tPrepareStateflow = 0;

double SchedulabilityAnalysis::tCalCSum = 0;
double SchedulabilityAnalysis::tCalLinearBounds = 0;
double SchedulabilityAnalysis::tCalTF0 = 0;
double SchedulabilityAnalysis::tCalTF1 = 0;
double SchedulabilityAnalysis::tCalTF2 = 0;

double SchedulabilityAnalysis::tCalNonPeriod = 0;
double SchedulabilityAnalysis::tCalPeriod = 0;
double SchedulabilityAnalysis::tCalDiffPeriod = 0;

double SchedulabilityAnalysis::tCalRequestExecutionMatrixWithoutPeriodicityProperty = 0;
double SchedulabilityAnalysis::tCalRequestExecutionMatrixWithPeriodicityProperty = 0;

double SchedulabilityAnalysis::tGenerateCriticalActionPairs = 0;
double SchedulabilityAnalysis::tGenerateRequestFunctions = 0;
double SchedulabilityAnalysis::tGenerateRequestFunctionAbstractTree = 0;

double SchedulabilityAnalysis::tExactStaticOffset = 0;
double SchedulabilityAnalysis::tRBFStaticOffset = 0;
double SchedulabilityAnalysis::tIBFStaticOffset = 0;
double SchedulabilityAnalysis::tRBFArbitraryOffset = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffset = 0;
double SchedulabilityAnalysis::tExactArbitraryOffsetBySimpleDigraph = 0; // Exact schedulability analysis on (simple) digraphs transformed with action instances
double SchedulabilityAnalysis::tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
double SchedulabilityAnalysis::tRBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
double SchedulabilityAnalysis::tExactArbitraryOffsetByPreciseDigraph = 0; // Exact schedulability analysis on (precise) digraphs transformed with action instances
double SchedulabilityAnalysis::tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
double SchedulabilityAnalysis::tRBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
double SchedulabilityAnalysis::tIBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
double SchedulabilityAnalysis::tLinearUpperRBF = 0;
double SchedulabilityAnalysis::tLinearUpperIBF = 0;
double SchedulabilityAnalysis::tLinearUpperCSUM = 0;

vector<double> SchedulabilityAnalysis::vec_tfRatio0;
vector<double> SchedulabilityAnalysis::vec_tfRatio1;

vector<int> SchedulabilityAnalysis::vec_nDigraphExact;
vector<int> SchedulabilityAnalysis::vec_nDigraphIBF;
vector<int> SchedulabilityAnalysis::vec_nDigraphRBF;

vector<int> SchedulabilityAnalysis::vec_nIrreducibleDigraphs;
vector<int> SchedulabilityAnalysis::vec_nDigraphNode;
vector<int> SchedulabilityAnalysis::vec_nDigraphEdge;

vector<double> SchedulabilityAnalysis::vec_tDigraphCalLinearFactor;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalLinearPeriod;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalLinearDefect;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalSum;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalLinearBounds;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalTF0;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalTF1;
vector<double> SchedulabilityAnalysis::vec_tDigraphCalTF2;

vector<double> SchedulabilityAnalysis::vec_tDigraphExact;
vector<double> SchedulabilityAnalysis::vec_tDigraphRBF;
vector<double> SchedulabilityAnalysis::vec_tDigraphIBF;

vector<int> SchedulabilityAnalysis::vec_nStateflows;
vector<double> SchedulabilityAnalysis::vec_totalUtilization;
vector<double> SchedulabilityAnalysis::vec_totalUtilization2;
vector<double> SchedulabilityAnalysis::vec_avgDegree;

vector<int> SchedulabilityAnalysis::vec_nTimeout;
vector<int> SchedulabilityAnalysis::vec_nExactStaticOffset;
vector<int> SchedulabilityAnalysis::vec_nRBFStaticOffset;
vector<int> SchedulabilityAnalysis::vec_nIBFStaticOffset;
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffset;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffset;
vector<int> SchedulabilityAnalysis::vec_nExactArbitraryOffsetBySimpleDigraph; // Exact schedulability analysis on (simple) digraphs transformed with action instances
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffsetBySimpleDigraphWithLinearization;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffsetBySimpleDigraphWithLinearization;
vector<int> SchedulabilityAnalysis::vec_nExactArbitraryOffsetByPreciseDigraph; // Exact schedulability analysis on (precise) digraphs transformed with action instances
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
vector<int> SchedulabilityAnalysis::vec_nRBFArbitraryOffsetByPreciseDigraphWithLinearization;
vector<int> SchedulabilityAnalysis::vec_nIBFArbitraryOffsetByPreciseDigraphWithLinearization;
vector<int> SchedulabilityAnalysis::vec_nLinearUpperRBF;
vector<int> SchedulabilityAnalysis::vec_nLinearUpperIBF;
vector<int> SchedulabilityAnalysis::vec_nLinearUpperCSUM;

vector<double> SchedulabilityAnalysis::vec_bpRatio0;
vector<double> SchedulabilityAnalysis::vec_bpRatio1;

vector<int> SchedulabilityAnalysis::vec_nIrreducibleStateflows;

vector<double> SchedulabilityAnalysis::vec_tPrepareStateflow;

vector<double> SchedulabilityAnalysis::vec_tCalCSum;
vector<double> SchedulabilityAnalysis::vec_tCalLinearBounds;
vector<double> SchedulabilityAnalysis::vec_tCalTF0;
vector<double> SchedulabilityAnalysis::vec_tCalTF1;
vector<double> SchedulabilityAnalysis::vec_tCalTF2;

vector<double> SchedulabilityAnalysis::vec_tCalNonPeriod;
vector<double> SchedulabilityAnalysis::vec_tCalPeriod;
vector<double> SchedulabilityAnalysis::vec_tCalDiffPeriod;

vector<double> SchedulabilityAnalysis::vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty;
vector<double> SchedulabilityAnalysis::vec_tCalRequestExecutionMatrixWithPeriodicityProperty;

vector<double> SchedulabilityAnalysis::vec_tGenerateCriticalActionPairs;
vector<double> SchedulabilityAnalysis::vec_tGenerateRequestFunctions;
vector<double> SchedulabilityAnalysis::vec_tGenerateRequestFunctionAbstractTree;

vector<double> SchedulabilityAnalysis::vec_tExactStaticOffset;
vector<double> SchedulabilityAnalysis::vec_tRBFStaticOffset;
vector<double> SchedulabilityAnalysis::vec_tIBFStaticOffset;
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffset;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffset;
vector<double> SchedulabilityAnalysis::vec_tExactArbitraryOffsetBySimpleDigraph; // Exact schedulability analysis on (simple) digraphs transformed with action instances
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph;
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffsetBySimpleDigraphWithLinearization;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffsetBySimpleDigraphWithLinearization;
vector<double> SchedulabilityAnalysis::vec_tExactArbitraryOffsetByPreciseDigraph; // Exact schedulability analysis on (precise) digraphs transformed with action instances
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph;
vector<double> SchedulabilityAnalysis::vec_tRBFArbitraryOffsetByPreciseDigraphWithLinearization;
vector<double> SchedulabilityAnalysis::vec_tIBFArbitraryOffsetByPreciseDigraphWithLinearization;
vector<double> SchedulabilityAnalysis::vec_tLinearUpperRBF;
vector<double> SchedulabilityAnalysis::vec_tLinearUpperIBF;
vector<double> SchedulabilityAnalysis::vec_tLinearUpperCSUM;

void SchedulabilityAnalysis::save_results() {
	vec_tfRatio0.push_back(tfRatio0);
	vec_tfRatio1.push_back(tfRatio1);

	vec_nDigraphExact.push_back(nDigraphExact);
	vec_nDigraphRBF.push_back(nDigraphRBF);
	vec_nDigraphIBF.push_back(nDigraphIBF);

	vec_nIrreducibleDigraphs.push_back(nIrreducibleDigraphs);
	vec_nDigraphNode.push_back(nDigraphNode);
	vec_nDigraphEdge.push_back(nDigraphEdge);

	vec_tDigraphCalLinearFactor.push_back(tDigraphCalLinearFactor);
	vec_tDigraphCalLinearPeriod.push_back(tDigraphCalLinearPeriod);
	vec_tDigraphCalLinearDefect.push_back(tDigraphCalLinearDefect);
	vec_tDigraphCalSum.push_back(tDigraphCalSum);
	vec_tDigraphCalLinearBounds.push_back(tDigraphCalLinearBounds);
	vec_tDigraphCalTF0.push_back(tDigraphCalTF0);
	vec_tDigraphCalTF1.push_back(tDigraphCalTF1);
	vec_tDigraphCalTF2.push_back(tDigraphCalTF2);

	vec_tDigraphExact.push_back(tDigraphExact);
	vec_tDigraphRBF.push_back(tDigraphRBF);
	vec_tDigraphIBF.push_back(tDigraphIBF);

	vec_nStateflows.push_back(nStateflows);
	vec_totalUtilization.push_back(totalUtilization);
	vec_totalUtilization2.push_back(totalUtilization2);
	vec_avgDegree.push_back(avgDegree);
	
	vec_nTimeout.push_back(nTimeout);
	vec_nExactStaticOffset.push_back(nExactStaticOffset);
	vec_nRBFStaticOffset.push_back(nRBFStaticOffset);
	vec_nIBFStaticOffset.push_back(nIBFStaticOffset);
	vec_nRBFArbitraryOffset.push_back(nRBFArbitraryOffset);
	vec_nIBFArbitraryOffset.push_back(nIBFArbitraryOffset);
	vec_nExactArbitraryOffsetBySimpleDigraph.push_back(nExactArbitraryOffsetBySimpleDigraph); // Exact schedulability analysis on (simple) digraphs transformed with action instances
	vec_nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph.push_back(nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	vec_nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph.push_back(nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	vec_nRBFArbitraryOffsetBySimpleDigraphWithLinearization.push_back(nRBFArbitraryOffsetBySimpleDigraphWithLinearization);
	vec_nIBFArbitraryOffsetBySimpleDigraphWithLinearization.push_back(nIBFArbitraryOffsetBySimpleDigraphWithLinearization);
	vec_nExactArbitraryOffsetByPreciseDigraph.push_back(nExactArbitraryOffsetByPreciseDigraph); // Exact schedulability analysis on (precise) digraphs transformed with action instances
	vec_nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph.push_back(nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	vec_nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph.push_back(nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	vec_nRBFArbitraryOffsetByPreciseDigraphWithLinearization.push_back(nRBFArbitraryOffsetByPreciseDigraphWithLinearization);
	vec_nIBFArbitraryOffsetByPreciseDigraphWithLinearization.push_back(nIBFArbitraryOffsetByPreciseDigraphWithLinearization);
	vec_nLinearUpperRBF.push_back(nLinearUpperRBF);
	vec_nLinearUpperIBF.push_back(nLinearUpperIBF);
	vec_nLinearUpperCSUM.push_back(nLinearUpperCSUM);

	vec_bpRatio0.push_back(bpRatio0);
	vec_bpRatio1.push_back(bpRatio1);

	vec_nIrreducibleStateflows.push_back(nIrreducibleStateflows);

	vec_tPrepareStateflow.push_back(tPrepareStateflow);

	vec_tCalCSum.push_back(tCalCSum);
	vec_tCalLinearBounds.push_back(tCalLinearBounds);
	vec_tCalTF0.push_back(tCalTF0);
	vec_tCalTF1.push_back(tCalTF1);
	vec_tCalTF2.push_back(tCalTF2);

	vec_tCalNonPeriod.push_back(tCalNonPeriod);
	vec_tCalPeriod.push_back(tCalPeriod);
	vec_tCalDiffPeriod.push_back(tCalDiffPeriod);

	vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty.push_back(tCalRequestExecutionMatrixWithoutPeriodicityProperty);
	vec_tCalRequestExecutionMatrixWithPeriodicityProperty.push_back(tCalRequestExecutionMatrixWithPeriodicityProperty);

	vec_tGenerateCriticalActionPairs.push_back(tGenerateCriticalActionPairs);
	vec_tGenerateRequestFunctions.push_back(tGenerateRequestFunctions);
	vec_tGenerateRequestFunctionAbstractTree.push_back(tGenerateRequestFunctionAbstractTree);
	
	vec_tExactStaticOffset.push_back(tExactStaticOffset);
	vec_tRBFStaticOffset.push_back(tRBFStaticOffset);
	vec_tIBFStaticOffset.push_back(tIBFStaticOffset);
	vec_tRBFArbitraryOffset.push_back(tRBFArbitraryOffset);
	vec_tIBFArbitraryOffset.push_back(tIBFArbitraryOffset);
	vec_tExactArbitraryOffsetBySimpleDigraph.push_back(tExactArbitraryOffsetBySimpleDigraph); // Exact schedulability analysis on (simple) digraphs transformed with action instances
	vec_tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph.push_back(tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	vec_tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph.push_back(tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	vec_tRBFArbitraryOffsetBySimpleDigraphWithLinearization.push_back(tRBFArbitraryOffsetBySimpleDigraphWithLinearization);
	vec_tIBFArbitraryOffsetBySimpleDigraphWithLinearization.push_back(tIBFArbitraryOffsetBySimpleDigraphWithLinearization);
	vec_tExactArbitraryOffsetByPreciseDigraph.push_back(tExactArbitraryOffsetByPreciseDigraph); // Exact schedulability analysis on (precise) digraphs transformed with action instances
	vec_tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph.push_back(tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	vec_tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph.push_back(tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	vec_tRBFArbitraryOffsetByPreciseDigraphWithLinearization.push_back(tRBFArbitraryOffsetByPreciseDigraphWithLinearization);
	vec_tIBFArbitraryOffsetByPreciseDigraphWithLinearization.push_back(tIBFArbitraryOffsetByPreciseDigraphWithLinearization);
	vec_tLinearUpperRBF.push_back(tLinearUpperRBF);
	vec_tLinearUpperIBF.push_back(tLinearUpperIBF);
	vec_tLinearUpperCSUM.push_back(tLinearUpperCSUM);
}

void SchedulabilityAnalysis::set_zero() {
	// reset all the parameters (=0)
	tfRatio0 = 0;
	tfRatio1 = 0;

	nDigraphExact = 0;
	nDigraphIBF = 0;
	nDigraphRBF = 0;

	nIrreducibleDigraphs =0;
	nDigraphNode = 0;
	nDigraphEdge = 0;

	tDigraphCalLinearFactor = 0;
	tDigraphCalLinearPeriod = 0;
	tDigraphCalLinearDefect = 0;
	tDigraphCalSum = 0;
	tDigraphCalLinearBounds = 0;
	tDigraphCalTF0 = 0;
	tDigraphCalTF1 = 0;
	tDigraphCalTF2 = 0;

	tDigraphExact = 0;
	tDigraphRBF = 0;
	tDigraphIBF = 0;

	nStateflows = 0;
	totalUtilization = 0;
	totalUtilization2 = 0;
	avgDegree = 0;

	nTimeout = 0;
	nExactStaticOffset = 0;
	nRBFStaticOffset = 0;
	nIBFStaticOffset = 0;
	nRBFArbitraryOffset = 0;
	nIBFArbitraryOffset = 0;
	nExactArbitraryOffsetBySimpleDigraph = 0; // Exact schedulability analysis on (simple) digraphs transformed with action instances
	nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
	nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
	nRBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
	nIBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
	nExactArbitraryOffsetByPreciseDigraph = 0; // Exact schedulability analysis on (precise) digraphs transformed with action instances
	nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
	nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
	nRBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
	nIBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
	nLinearUpperRBF = 0;
	nLinearUpperIBF = 0;
	nLinearUpperCSUM = 0;

	bpRatio0 = 0;
	bpRatio1 = 0;

	nIrreducibleStateflows = 0;

	tPrepareStateflow = 0;

	tCalCSum = 0;
	tCalLinearBounds = 0;
	tCalTF0 = 0;
	tCalTF1 = 0;
	tCalTF2 = 0;

	tCalNonPeriod = 0;
	tCalPeriod = 0;
	tCalDiffPeriod = 0;

	tCalRequestExecutionMatrixWithoutPeriodicityProperty = 0;
	tCalRequestExecutionMatrixWithPeriodicityProperty = 0;

	tGenerateCriticalActionPairs = 0;
	tGenerateRequestFunctions = 0;
	tGenerateRequestFunctionAbstractTree = 0;

	tExactStaticOffset = 0;
	tRBFStaticOffset = 0;
	tIBFStaticOffset = 0;
	tRBFArbitraryOffset = 0;
	tIBFArbitraryOffset = 0;
	tExactArbitraryOffsetBySimpleDigraph = 0; // Exact schedulability analysis on (simple) digraphs transformed with action instances
	tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
	tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph = 0;
	tRBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
	tIBFArbitraryOffsetBySimpleDigraphWithLinearization = 0;
	tExactArbitraryOffsetByPreciseDigraph = 0; // Exact schedulability analysis on (precise) digraphs transformed with action instances
	tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
	tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph = 0;
	tRBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
	tIBFArbitraryOffsetByPreciseDigraphWithLinearization = 0;
	tLinearUpperRBF = 0;
	tLinearUpperIBF = 0;
	tLinearUpperCSUM = 0;
}

void SchedulabilityAnalysis::clear_all() {
	set_zero();

	vec_tfRatio0.clear();
	vec_tfRatio1.clear();

	vec_nDigraphExact.clear();
	vec_nDigraphIBF.clear();
	vec_nDigraphRBF.clear();

	vec_nIrreducibleDigraphs.clear();
	vec_nDigraphNode.clear();
	vec_nDigraphEdge.clear();

	vec_tDigraphCalLinearFactor.clear();
	vec_tDigraphCalLinearPeriod.clear();
	vec_tDigraphCalLinearDefect.clear();
	vec_tDigraphCalSum.clear();
	vec_tDigraphCalLinearBounds.clear();
	vec_tDigraphCalTF0.clear();
	vec_tDigraphCalTF1.clear();
	vec_tDigraphCalTF2.clear();

	vec_tDigraphExact.clear();
	vec_tDigraphRBF.clear();
	vec_tDigraphIBF.clear();

	vec_nStateflows.clear();
	vec_totalUtilization.clear();
	vec_totalUtilization2.clear();
	vec_avgDegree.clear();

	vec_nTimeout.clear();
	vec_nExactStaticOffset.clear();
	vec_nRBFStaticOffset.clear();
	vec_nIBFStaticOffset.clear();
	vec_nRBFArbitraryOffset.clear();
	vec_nIBFArbitraryOffset.clear();
	vec_nExactArbitraryOffsetBySimpleDigraph.clear(); // Exact schedulability analysis on (simple) digraphs transformed with action instances
	vec_nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph.clear();
	vec_nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph.clear();
	vec_nRBFArbitraryOffsetBySimpleDigraphWithLinearization.clear();
	vec_nIBFArbitraryOffsetBySimpleDigraphWithLinearization.clear();
	vec_nExactArbitraryOffsetByPreciseDigraph.clear(); // Exact schedulability analysis on (precise) digraphs transformed with action instances
	vec_nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph.clear();
	vec_nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph.clear();
	vec_nRBFArbitraryOffsetByPreciseDigraphWithLinearization.clear();
	vec_nIBFArbitraryOffsetByPreciseDigraphWithLinearization.clear();
	vec_nLinearUpperRBF.clear();
	vec_nLinearUpperIBF.clear();
	vec_nLinearUpperCSUM.clear();

	vec_bpRatio0.clear();
	vec_bpRatio1.clear();
	vec_nIrreducibleStateflows.clear();

	vec_tPrepareStateflow.clear();

	vec_tCalCSum.clear();
	vec_tCalLinearBounds.clear();
	vec_tCalTF0.clear();
	vec_tCalTF1.clear();
	vec_tCalTF2.clear();
	vec_tCalNonPeriod.clear();
	vec_tCalPeriod.clear();
	vec_tCalDiffPeriod.clear();
	vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty.clear();
	vec_tCalRequestExecutionMatrixWithPeriodicityProperty.clear();
	vec_tGenerateCriticalActionPairs.clear();
	vec_tGenerateRequestFunctions.clear();
	vec_tGenerateRequestFunctionAbstractTree.clear();

	vec_tExactStaticOffset.clear();
	vec_tRBFStaticOffset.clear();
	vec_tIBFStaticOffset.clear();
	vec_tRBFArbitraryOffset.clear();
	vec_tIBFArbitraryOffset.clear();
	vec_tExactArbitraryOffsetBySimpleDigraph.clear(); // Exact schedulability analysis on (simple) digraphs transformed with action instances
	vec_tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph.clear();
	vec_tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph.clear();
	vec_tRBFArbitraryOffsetBySimpleDigraphWithLinearization.clear();
	vec_tIBFArbitraryOffsetBySimpleDigraphWithLinearization.clear();
	vec_tExactArbitraryOffsetByPreciseDigraph.clear(); // Exact schedulability analysis on (precise) digraphs transformed with action instances
	vec_tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph.clear();
	vec_tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph.clear();
	vec_tRBFArbitraryOffsetByPreciseDigraphWithLinearization.clear();
	vec_tIBFArbitraryOffsetByPreciseDigraphWithLinearization.clear();
	vec_tLinearUpperRBF.clear();
	vec_tLinearUpperIBF.clear();
	vec_tLinearUpperCSUM.clear();
}

void SchedulabilityAnalysis::reset() {
	// save the results into the vectors
	save_results();
	set_zero();
}

void SchedulabilityAnalysis::output_one_vector(ofstream& fout, string sVec, vector<double> vec) {
	fout<<sVec<<"=[";
	typedef vector<double>::iterator Iter;
	for (Iter it = vec.begin(); it != vec.end(); it++) {
		fout<<*it;
		if (it != --vec.end()) fout<<",";
	}
	fout<<"];"<<endl;
}

void SchedulabilityAnalysis::output_one_vector(ofstream& fout, string sVec, vector<int> vec) {
	fout<<sVec<<"=[";
	typedef vector<int>::iterator Iter;
	for (Iter it = vec.begin(); it != vec.end(); it++) {
		fout<<*it;
		if (it != --vec.end()) fout<<",";
	}
	fout<<"];"<<endl;
}

void SchedulabilityAnalysis::output_vectors(ofstream& fout)  {
	fout<<endl;
	fout<<"Output the informations of Digraphs"<<endl;
	fout<<endl;

	output_one_vector(fout,"tfRatio0",vec_tfRatio0);
	output_one_vector(fout,"tfRatio1",vec_tfRatio1);

	output_one_vector(fout,"nDigraphExact",vec_nDigraphExact);
	output_one_vector(fout,"nDigraphRBF",vec_nDigraphRBF);
	output_one_vector(fout,"nDigraphIBF",vec_nDigraphIBF);

	output_one_vector(fout,"nIrreducibleDigraphs",vec_nIrreducibleDigraphs);
	output_one_vector(fout,"nDigraphNode",vec_nDigraphNode);
	output_one_vector(fout,"nDigraphEdge",vec_nDigraphEdge);

	output_one_vector(fout,"tDigraphCalLinearFactor",vec_tDigraphCalLinearFactor);
	output_one_vector(fout,"tDigraphCalLinearPeriod",vec_tDigraphCalLinearPeriod);
	output_one_vector(fout,"tDigraphCalLinearDefect",vec_tDigraphCalLinearDefect);
	output_one_vector(fout,"tDigraphCalSum",vec_tDigraphCalSum);
	output_one_vector(fout,"tDigraphCalLinearBounds",vec_tDigraphCalLinearBounds);
	output_one_vector(fout,"tDigraphCalTF0",vec_tDigraphCalTF0);
	output_one_vector(fout,"tDigraphCalTF1",vec_tDigraphCalTF1);
	output_one_vector(fout,"tDigraphCalTF2",vec_tDigraphCalTF2);

	output_one_vector(fout,"tDigraphExact",vec_tDigraphExact);
	output_one_vector(fout,"tDigraphRBF",vec_tDigraphRBF);
	output_one_vector(fout,"tDigraphIBF",vec_tDigraphIBF);

	fout<<endl;
	fout<<"Output the informations of Stateflows"<<endl;
	fout<<endl;

	output_one_vector(fout,"nStateflows",vec_nStateflows);
	output_one_vector(fout,"totalUtil",vec_totalUtilization);
	output_one_vector(fout,"totalUtil2",vec_totalUtilization2);
	output_one_vector(fout,"avgDegree",vec_avgDegree);

	output_one_vector(fout,"nTimeout",vec_nTimeout);
	output_one_vector(fout,"nExactStaticOffset", vec_nExactStaticOffset);
	output_one_vector(fout,"nRBFStaticOffset",vec_nRBFStaticOffset);
	output_one_vector(fout,"nIBFStaticOffset",vec_nIBFStaticOffset);
	output_one_vector(fout,"nRBFArbitraryOffset",vec_nRBFArbitraryOffset);
	output_one_vector(fout,"nIBFArbitraryOffset",vec_nIBFArbitraryOffset);

	output_one_vector(fout,"nExactArbitraryOffsetBySimpleDigraph",vec_nExactArbitraryOffsetBySimpleDigraph); 
	output_one_vector(fout,"nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph", vec_nRBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	output_one_vector(fout,"nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph", vec_nIBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	output_one_vector(fout,"nRBFArbitraryOffsetBySimpleDigraphWithLinearization", vec_nRBFArbitraryOffsetBySimpleDigraphWithLinearization);
	output_one_vector(fout,"nIBFArbitraryOffsetBySimpleDigraphWithLinearization", vec_nIBFArbitraryOffsetBySimpleDigraphWithLinearization);
	
	output_one_vector(fout,"nExactArbitraryOffsetByPreciseDigraph",vec_nExactArbitraryOffsetByPreciseDigraph); 
	output_one_vector(fout,"nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph", vec_nRBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	output_one_vector(fout,"nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph", vec_nIBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	output_one_vector(fout,"nRBFArbitraryOffsetByPreciseDigraphWithLinearization", vec_nRBFArbitraryOffsetByPreciseDigraphWithLinearization);
	output_one_vector(fout,"nIBFArbitraryOffsetByPreciseDigraphWithLinearization", vec_nIBFArbitraryOffsetByPreciseDigraphWithLinearization);

	output_one_vector(fout,"nLinearUpperRBF",vec_nLinearUpperRBF);
	output_one_vector(fout,"nLinearUpperIBF",vec_nLinearUpperIBF);
	output_one_vector(fout,"nLinearUpperCSUM", vec_nLinearUpperCSUM);

	output_one_vector(fout,"bpRatio0",vec_bpRatio0);
	output_one_vector(fout,"bpRatio1",vec_bpRatio1);

	output_one_vector(fout,"nIrreducibleStateflows",vec_nIrreducibleStateflows);

	output_one_vector(fout,"tPrepareStateflow",vec_tPrepareStateflow);

	output_one_vector(fout,"tCalCSum",vec_tCalCSum);
	output_one_vector(fout,"tCalLinearBounds",vec_tCalLinearBounds);
	output_one_vector(fout,"tCalTF0",vec_tCalTF0);
	output_one_vector(fout,"tCalTF1",vec_tCalTF1);
	output_one_vector(fout,"tCalTF2",vec_tCalTF2);

	output_one_vector(fout,"tCalNonPeirod",vec_tCalNonPeriod);
	output_one_vector(fout,"tCalPeriod",vec_tCalPeriod);
	output_one_vector(fout,"tCalDiffPeriod",vec_tCalDiffPeriod);

	output_one_vector(fout,"tCalRequestExecutionMatrixWithoutPeriodicityProperty",vec_tCalRequestExecutionMatrixWithoutPeriodicityProperty);
	output_one_vector(fout,"tCalRequestExecutionMatrixWithPeriodicityProperty",vec_tCalRequestExecutionMatrixWithPeriodicityProperty);

	output_one_vector(fout,"tGenerateCriticalPairs",vec_tGenerateCriticalActionPairs);
	output_one_vector(fout,"tGenerateRequestFunctions", vec_tGenerateRequestFunctions);
	output_one_vector(fout,"tGenerateRequestFunctionAbstractTree", vec_tGenerateRequestFunctionAbstractTree);
	
	output_one_vector(fout,"tExactStaticOffset", vec_tExactStaticOffset);
	output_one_vector(fout,"tRBFStaticOffset",vec_tRBFStaticOffset);
	output_one_vector(fout,"tIBFStaticOffset",vec_tIBFStaticOffset);
	output_one_vector(fout,"tRBFArbitraryOffset",vec_tRBFArbitraryOffset);
	output_one_vector(fout,"tIBFArbitraryOffset",vec_tIBFArbitraryOffset);

	output_one_vector(fout,"tExactArbitraryOffsetBySimpleDigraph",vec_tExactArbitraryOffsetBySimpleDigraph); 
	output_one_vector(fout,"tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph", vec_tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	output_one_vector(fout,"tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph", vec_tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph);
	output_one_vector(fout,"tRBFArbitraryOffsetBySimpleDigraphWithLinearization", vec_tRBFArbitraryOffsetBySimpleDigraphWithLinearization);
	output_one_vector(fout,"tIBFArbitraryOffsetBySimpleDigraphWithLinearization", vec_tIBFArbitraryOffsetBySimpleDigraphWithLinearization);
	
	output_one_vector(fout,"tExactArbitraryOffsetByPreciseDigraph",vec_tExactArbitraryOffsetByPreciseDigraph); 
	output_one_vector(fout,"tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph", vec_tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	output_one_vector(fout,"tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph", vec_tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph);
	output_one_vector(fout,"tRBFArbitraryOffsetByPreciseDigraphWithLinearization", vec_tRBFArbitraryOffsetByPreciseDigraphWithLinearization);
	output_one_vector(fout,"tIBFArbitraryOffsetByPreciseDigraphWithLinearization", vec_tIBFArbitraryOffsetByPreciseDigraphWithLinearization);

	output_one_vector(fout,"tLinearUpperRBF",vec_tLinearUpperRBF);
	output_one_vector(fout,"tLinearUpperIBF",vec_tLinearUpperIBF);
	output_one_vector(fout,"tLinearUpperCSUM", vec_tLinearUpperCSUM);
}

void SchedulabilityAnalysis::prepare_all_digraphs(Digraph** digraphs, int n, int choice) {
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
		Timer timer;
		timer.start();
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		digraph->calculate_linear_factor();
		timer.end();
		tDigraphCalLinearFactor += timer.getTime();

		if (choice == 0) {
			timer.start();
			digraph->calculate_csum();
			timer.end();
			tDigraphCalSum += timer.getTime();

			timer.start();
			digraph->calculate_tf0(digraphs,i);
			timer.end();
			tDigraphCalTF0 += timer.getTime();
		}

		// Linear upper bounds should be finished in Digraph::prepare_digraphs
		if (choice == 1 || choice == 2) {
			timer.start();
			digraph->calculate_linear_upper_bounds();
			timer.end();
			tDigraphCalLinearBounds += timer.getTime();
		

			timer.start();
			digraph->calculate_tf1(digraphs,i);
			timer.end();
			tDigraphCalTF1 += timer.getTime();

			timer.start();
			digraph->calculate_tf2(digraphs,i);
			timer.end();
			tDigraphCalTF2 += timer.getTime();
		
			tfRatio0 += digraph->tf1/digraph->tf0;
			tfRatio1 += digraph->tf2/digraph->tf0;
		}
	}
}

void SchedulabilityAnalysis::prepare_all_digraphs2(Digraph** digraphs, int n) {
	Timer timer;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();

		timer.start();
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		digraph->calculate_linear_factor();
		timer.end();
		tDigraphCalLinearFactor += timer.getTime();

		timer.start();
		digraph->calculate_csum();
		timer.end();
		tDigraphCalSum += timer.getTime();

		timer.start();
		digraph->calculate_tf0(digraphs,i);
		timer.end();
		tDigraphCalTF0 += timer.getTime();

		timer.start();
		digraph->calculate_linear_upper_bounds();
		timer.end();
		tDigraphCalLinearBounds += timer.getTime();

		timer.start();
		digraph->calculate_tf1(digraphs,i);
		timer.end();
		tDigraphCalTF1 += timer.getTime();

		timer.start();
		digraph->calculate_tf2(digraphs,i);
		timer.end();
		tDigraphCalTF2 += timer.getTime();
		
		tfRatio0 += digraph->tf1/digraph->tf0;
		tfRatio1 += digraph->tf2/digraph->tf0;
	}
}

bool SchedulabilityAnalysis::isAbstract(vector<DigraphRFNode*> arfs) {
	for (vector<DigraphRFNode*>::iterator iter = arfs.begin(); iter != arfs.end(); iter++) {
		DigraphRFNode* rfnode = *iter;
		if (rfnode->isAbstract) return true;
	}
	return false;
}

void SchedulabilityAnalysis::generate_critical_vertices(Digraph** digraphs, int n) {
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		// Determine intra-digraph dominances
		set<Node*> Dominated;

		for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
			Node* node = *iter;
			if (Dominated.empty()) { 
				Dominated.insert(node);
				continue;
			}

			bool dominated = false;
			bool dominating = false;

			set<Node*> remove;

			for (set<Node*>::iterator sIter = Dominated.begin(); sIter != Dominated.end(); sIter++) {
				Node* cnode = *sIter;
				if (node->wcet >= cnode->wcet && node->deadline <= cnode->deadline) {
					dominating = true;
					remove.insert(cnode);
					continue;
				}

				if (node->wcet <= cnode->wcet && node->deadline >= cnode->deadline) {
					dominated = true;
				}
			}

			if (dominating) {
				for (set<Node*>::iterator sIter = remove.begin(); sIter != remove.end(); sIter++)
					Dominated.erase(*sIter);
			}

			if (!dominated)
				Dominated.insert(node);
		}

		for (set<Node*>::iterator iter = Dominated.begin(); iter != Dominated.end(); iter++) {
			digraph->cnode_vec.push_back(*iter);
		}

		// Determine inter-digraph dominances
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			set<Node*> Dominated2;
			vector<Node*> cnode_vec;
			for (vector<Node*>::iterator iiter = digraph->cnode_vec.begin(); iiter != digraph->cnode_vec.end(); iiter++) {
				Node* inode = *iiter;
				for (vector<Node*>::iterator jiter = digraphj->cnode_vec.begin(); jiter != digraphj->cnode_vec.end(); jiter++) {
					Node* jnode = *jiter;
					if (Dominated2.find(jnode) != Dominated2.end()) continue;
					if (inode->wcet >= jnode->wcet && inode->deadline <= jnode->deadline) Dominated2.insert(jnode);
				}
			}

			for (vector<Node*>::iterator iter = digraphj->cnode_vec.begin(); iter != digraphj->cnode_vec.end(); iter++) {
				Node* node = *iter;
				if (Dominated2.find(node) == Dominated2.end()) cnode_vec.push_back(node);
			}

			digraphj->cnode_vec = cnode_vec;
		}
	}
}

void SchedulabilityAnalysis::generate_critical_vertices_for_reachability_digraphs(Digraph** digraphs, int n) {
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		// Determine intra-digraph dominances
		set<Edge*> Dominated;

		for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
			Edge* edge = *iter;
			if (Dominated.empty()) { 
				Dominated.insert(edge);
				continue;
			}

			bool dominated = false;
			bool dominating = false;

			set<Edge*> remove;

			for (set<Edge*>::iterator sIter = Dominated.begin(); sIter != Dominated.end(); sIter++) {
				Edge* cedge = *sIter;
				if (edge->weight >= cedge->weight && edge->deadline <= cedge->deadline) {
					dominating = true;
					remove.insert(cedge);
					continue;
				}

				if (edge->weight <= cedge->weight && edge->deadline >= cedge->deadline) {
					dominated = true;
				}
			}

			if (dominating) {
				for (set<Edge*>::iterator sIter = remove.begin(); sIter != remove.end(); sIter++)
					Dominated.erase(*sIter);
			}

			if (!dominated)
				Dominated.insert(edge);
		}

		for (set<Edge*>::iterator iter = Dominated.begin(); iter != Dominated.end(); iter++) {
			digraph->cedge_vec.push_back(*iter);
		}

		// Determine inter-digraph dominances
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			set<Edge*> Dominated2;
			vector<Edge*> cedge_vec;
			for (vector<Edge*>::iterator iiter = digraph->cedge_vec.begin(); iiter != digraph->cedge_vec.end(); iiter++) {
				Edge* iedge = *iiter;
				for (vector<Edge*>::iterator jiter = digraphj->cedge_vec.begin(); jiter != digraphj->cedge_vec.end(); jiter++) {
					Edge* jedge = *jiter;
					if (Dominated2.find(jedge) != Dominated2.end()) continue;
					if (iedge->weight >= jedge->weight && iedge->deadline <= jedge->deadline) Dominated2.insert(jedge);
				}
			}

			for (vector<Edge*>::iterator iter = digraphj->cedge_vec.begin(); iter != digraphj->cedge_vec.end(); iter++) {
				Edge* edge = *iter;
				if (Dominated2.find(edge) == Dominated2.end()) cedge_vec.push_back(edge);
			}

			digraphj->cedge_vec = cedge_vec;
		}
	}
}

double SchedulabilityAnalysis::calculate_digraph_exact_response_time(Digraph** digraphs, int i, int wcet) {
	vector<DigraphRFNode*> combination;
	map<double,vector<DigraphRFNode*>> store;

	for (int j=0; j<i; j++) {
		Digraph* digraphj = digraphs[j];
		combination.push_back(digraphj->root);
	}

	double rt = calculate_digraph_exact_response_time(combination,wcet);
	store[rt] = combination;

	while (isAbstract((--store.end())->second)) {
		vector<DigraphRFNode*> combination1;
		vector<DigraphRFNode*> combination2;

		bool found = false;
		for (vector<DigraphRFNode*>::iterator iter = combination.begin(); iter != combination.end(); iter++) {
			DigraphRFNode* rfnode = *iter;
			if (rfnode->isAbstract && !found) {
				found = true;
				combination1.push_back(rfnode->lnode);
				combination2.push_back(rfnode->rnode);
				continue;
			}
			combination1.push_back(rfnode);
			combination2.push_back(rfnode);
		}

		double rt1 = calculate_digraph_exact_response_time(combination1,wcet);
		double rt2 = calculate_digraph_exact_response_time(combination2,wcet);

		// remove the last element
		store.erase(--store.end());
		store[rt1] = combination1;
		store[rt2] = combination2;
	}

	return (--store.end())->first;
}

double SchedulabilityAnalysis::calculate_digraph_exact_response_time(vector<DigraphRFNode*> combination, int wcet) {
	int gcd = 0;
	for (vector<DigraphRFNode*>::iterator iter = combination.begin(); iter != combination.end(); iter++) {
		gcd = Utility::math_gcd(gcd, (*iter)->gcd);
	}

	int ub = (*combination.begin())->ub;

	// check the schedulability condition
	for (int t=0; t<=ub; t+= gcd) {
		double sum = wcet;
		for (vector<DigraphRFNode*>::iterator iter = combination.begin(); iter != combination.end(); iter++) {
			DigraphRFNode* arf = *iter;
			sum += arf->get_mrf(t);
		}
		
		if (sum <= t) return sum;

		t = floor(1.0*sum/gcd)*gcd;
	}

	return 0;
}

bool SchedulabilityAnalysis::digraph_exact_analysis_timeout(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();

	if (false) cout << "Do analysis ..." << endl;
	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_exact_analysis(digraphs, i);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphExact += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_exact_analysis(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();

	//prepare_all_digraphs(digraphs,n,choice);

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
	}

	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	// prepare request functions
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i]; 
		if (false) cout << "Generate critical request functions..." << endl;
		digraph->generate_critical_request_function(tf_max,false);
		if (false) cout << "Generate abstract request functions tree..." << endl;
		digraph->generate_abstract_request_function_tree();
	}

	if (false) cout << "Do analysis ..." << endl;
	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_exact_analysis(digraphs, i);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphExact += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_exact_analysis(Digraph** digraphs, int i) {
	Digraph* digraphi = digraphs[i];

	vector<DigraphRFNode*> combination;
	for (int j=0; j<i; j++) {
		Digraph* digraphj = digraphs[j];
		combination.push_back(digraphj->root);
	}

	for(vector<Node*>::iterator iter = digraphi->cnode_vec.begin(); iter != digraphi->cnode_vec.end(); iter++) {
		Node* cnode = *iter;
		bool schedulable = digraph_exact_analysis(combination, cnode->wcet, cnode->deadline);
		if (!schedulable) return false;
	}
	
	return true;
}

bool SchedulabilityAnalysis::digraph_exact_analysis(vector<DigraphRFNode*> combination, int wcet, int deadline) {
	if (combination.empty()) {
		if (wcet <= deadline) return true;
		else return false;
	}
	
	int gcd = 0;
	for (vector<DigraphRFNode*>::iterator iter = combination.begin(); iter != combination.end(); iter++) {
		gcd = Utility::math_gcd(gcd, (*iter)->gcd);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+= gcd) {
		double sum = wcet;
		for (vector<DigraphRFNode*>::iterator iter = combination.begin(); iter != combination.end(); iter++) {
			DigraphRFNode* arf = *iter;
			sum += arf->get_mrf(t);
		}
		
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}

	if (!isAbstract(combination)) return false;

	vector<DigraphRFNode*> vec_rfnode1;
	vector<DigraphRFNode*> vec_rfnode2;

	bool found = false;
	for (vector<DigraphRFNode*>::iterator iter = combination.begin(); iter != combination.end(); iter++) {
		DigraphRFNode* rfnode = *iter;
		if (rfnode->isAbstract && !found) {
			found = true;
			vec_rfnode1.push_back(rfnode->lnode);
			vec_rfnode2.push_back(rfnode->rnode);
			continue;
		}
		vec_rfnode1.push_back(rfnode);
		vec_rfnode2.push_back(rfnode);
	}

	bool schedulable1 = digraph_exact_analysis(vec_rfnode1,wcet, deadline);
	bool schedulable2 = digraph_exact_analysis(vec_rfnode2,wcet, deadline);
		
	return schedulable1 && schedulable2;
}

/*
bool SchedulabilityAnalysis::digraph_exact_analysis(vector<DigraphRFNode*> roots, int wcet, int deadline) {
	vector<vector<DigraphRFNode*>> combinations;
	combinations.push_back(roots);
	
	return digraph_exact_analysis(combinations, wcet, deadline);
}

bool SchedulabilityAnalysis::digraph_exact_analysis(vector<vector<DigraphRFNode*>> combinations, int wcet, int deadline) {
	for (vector<vector<DigraphRFNode*>>::iterator iter = combinations.begin(); iter != combinations.end(); iter++) {
		vector<DigraphRFNode*> combination = *iter;
		bool schedulable = digraph_exact_analysis(combination,wcet, deadline);
		if (!schedulable) return false;
	}
	return true;
}
*/

bool SchedulabilityAnalysis::digraph_rbf_analysis(Digraph** digraphs, int n, int choice, bool linearization) {
	Timer timer;
	timer.start();

	if (linearization)
		prepare_all_digraphs(digraphs,n,choice);
	else {
		for (int i=0; i<n; i++) {
			Digraph* digraph = digraphs[i];
			digraph->calculate_period_gcd();
			digraph->calculate_all_gcd();
		}
	}

	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;

		if (linearization) {
			//digraph->prepare_rbf_calculation(false);
			//int tf = ceil(1.0*tf_max/digraph->pGCD)+10;
			//digraph->unit_digraph->calculate_exec_request_matrix_power(tf);
			digraph->prepare_rbf_calculation_without_periodicity(false);
		}
		else {
			int tf = ceil(1.0*tf_max/digraph->pGCD)*digraph->pGCD;
			digraph->calculate_rbf_without_periodicity(tf,digraph->rbf_map2);
			//digraph->calculate_rbf_without_periodicity_fast(tf,digraph->rbf_map2_fast);
		}
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_rbf_analysis(digraphs,i,linearization);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphRBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis(Digraph** digraphs, int i, bool linearization) {
	Digraph* digraph = digraphs[i];
	for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
		Node* node = *iter;
		bool schedulable = digraph_rbf_analysis(digraphs, i, linearization, node->wcet, node->deadline);
		if (!schedulable) {
			cout << "Digraph-" << i << ": (wcet,deadline)=<" <<node->wcet<< "," << node->deadline << ">" << endl; 	
			return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis(Digraph** digraphs, int i, bool linearization, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->pGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			if (linearization) 
				sum += digraphj->rbf(t);
			else
				sum += digraphj->rbf2(t);
				//sum += digraphj->rbf2_fast(t);
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return false;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis2(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();

	prepare_all_digraphs(digraphs,n,choice);

	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;

		digraph->prepare_rbf_calculation(false);

		if (digraph->unit_digraph->ldef == -1) {
			digraph->calculate_rbf_without_periodicity(tf_max, digraph->rbf_map);
		} else {
			digraph->calculate_rbf_without_periodicity((digraph->unit_digraph->ldef+digraph->unit_digraph->lper)*digraph->pGCD, digraph->rbf_map);
			int lper = digraph->unit_digraph->lper*digraph->pGCD;
			for (int t=(digraph->unit_digraph->ldef+digraph->unit_digraph->lper)*digraph->pGCD+1; t <= tf_max; t++) {
				digraph->rbf_map[t] = digraph->rbf_map[t-lper] + digraph->linear_factor*lper;  
			}
		}
			
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_rbf_analysis2(digraphs,i);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphRBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis2(Digraph** digraphs, int i) {
	Digraph* digraph = digraphs[i];
	for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
		Node* node = *iter;
		bool schedulable = digraph_rbf_analysis2(digraphs, i, node->wcet, node->deadline);
		if (!schedulable) {
			cout << "Digraph-" << i << ": (wcet,deadline)=<" <<node->wcet<< "," << node->deadline << ">" << endl; 	
			return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis2(Digraph** digraphs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->pGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			map<int,int>::const_iterator mIter = digraphj->rbf_map.lower_bound(t);
			sum += mIter->second;
				//sum += digraphj->rbf2_fast(t);
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return false;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice,bool periodicityProperty, bool deadlineProperty) {
	Timer timer;
	timer.start();

	if (deadlineProperty) {
		prepare_all_digraphs(digraphs,n,2);

		for (int i=0; i<n; i++) {
			Digraph* digraph = digraphs[i];

			//cout << "digraph-"<< i << "=>tf=" << digraph->tf2 << endl;

			if (periodicityProperty && digraph->strongly_connected) {
				digraph->calculate_linear_factor();
				digraph->unit_digraph = new UnitDigraph(digraph);
				digraph->unit_digraph->prepare3(false);
			}
		}
	}
	else {
		for (int i=0; i<n; i++) {
			Digraph* digraph = digraphs[i];
			digraph->calculate_period_gcd();
			digraph->calculate_all_gcd();

			if (periodicityProperty) {
				digraph->generate_strongly_connected_components();
				digraph->check_strongly_connected();

				if (digraph->strongly_connected) {
					digraph->calculate_linear_factor();
					digraph->unit_digraph = new UnitDigraph(digraph);
					digraph->unit_digraph->prepare3(false);
				}
			}
		}
	}
	
	
	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		//nDigraphNode += digraph->node_vec.size();
		//nDigraphEdge += digraph->edge_vec.size();

		if (deadlineProperty)
			tf_max = max(tf_max,digraph->tf2);
		else {
			for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
				tf_max = max(tf_max,(*iter)->deadline);
			}
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;
		int tf = ceil(1.0*tf_max/digraph->pGCD)*digraph->pGCD;

		if (periodicityProperty && digraph->strongly_connected) {
			nIrreducibleDigraphs ++;
			digraph->calculate_rbf_with_periodicity_DP(tf,digraph->rbf_map2_DP);
		}
		else {
			digraph->calculate_rbf_without_periodicity_DP(tf,digraph->rbf_map2_DP);	
		}

		if (false) {
			cout << "==============" << "Digraph-" << i << "==============" << endl;
			for (map<int,int>::iterator mIter = digraph->rbf_map2_DP.begin(); mIter != digraph->rbf_map2_DP.end(); mIter++) 
				cout << "Digraph-" << i << ".rbf(" << mIter->first << ")=" <<mIter->second << endl;
		}
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_rbf_analysis_dynamic_programming(digraphs,i);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphRBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i) {
	Digraph* digraph = digraphs[i];
	for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
		Node* node = *iter;
		/*
		if (i==17) {
			cout << "Critical Vertex=(" << node->wcet << "," << node->deadline << endl;  
		}
		*/
		bool schedulable = digraph_rbf_analysis_dynamic_programming(digraphs, i, node->wcet, node->deadline);
		if (!schedulable) {
			cout << "Digraph-" << i << ": (wcet,deadline)=<" <<node->wcet<< "," << node->deadline << ">" << endl; 	
			return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->pGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->rbf2(t,digraphj->rbf_map2_DP);
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}

	return false;
}

bool SchedulabilityAnalysis::reachability_digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();

	
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
	}
	
	generate_critical_vertices_for_reachability_digraphs(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		//nDigraphNode += digraph->node_vec.size();
		//nDigraphEdge += digraph->edge_vec.size();

		for (vector<Edge*>::iterator iter = digraph->cedge_vec.begin(); iter != digraph->cedge_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;
		
		digraph->calculate_rbf_without_periodicity_on_reachability_digraph_DP(tf_max,digraph->rbf_map2_DP);	
		
		if (false) {
			cout << "==============" << "Digraph-" << i << "==============" << endl;
			for (map<int,int>::iterator mIter = digraph->rbf_map2_DP.begin(); mIter != digraph->rbf_map2_DP.end(); mIter++) 
				cout << "Digraph-" << i << ".rbf(" << mIter->first << ")=" <<mIter->second << endl;
		}
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::reachability_digraph_rbf_analysis_dynamic_programming(digraphs,i);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphRBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::reachability_digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i) {
	Digraph* digraph = digraphs[i];
	for (vector<Edge*>::iterator iter = digraph->cedge_vec.begin(); iter != digraph->cedge_vec.end(); iter++) {
		Edge* edge = *iter;
		/*
		if (i==17) {
			cout << "Critical Vertex=(" << edge->weight << "," << edge->deadline << endl;  
		}
		*/
		bool schedulable = reachability_digraph_rbf_analysis_dynamic_programming(digraphs, i, edge->weight, edge->deadline);
		if (!schedulable) {
			cout << "Digraph-" << i << ": (wcet,deadline)=<" <<edge->weight<< "," << edge->deadline << ">" << endl; 	
			return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::reachability_digraph_rbf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->pGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->rbf2(t,digraphj->rbf_map2_DP);
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}

	return false;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis(Digraph** digraphs, int n, int choice, bool linearization) {
	Timer timer;
	timer.start();

	if (linearization)
		prepare_all_digraphs(digraphs,n,choice);
	else {
		for (int i=0; i<n; i++) {
			Digraph* digraph = digraphs[i];
			digraph->calculate_period_gcd();
			digraph->calculate_all_gcd();
		}
	}

	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;

		if (linearization) 
			digraph->prepare_ibf_calculation2(false);
		else {
			int tf = ceil(1.0*tf_max/digraph->aGCD)*digraph->aGCD;
			//digraph->calculate_ibf_without_periodicity(tf,digraph->ibf_map2);
			digraph->calculate_ibf_without_periodicity_fast(tf,digraph->ibf_map2_fast);

			/*
			if (true) {
				for (int t=0; t<=tf; t+= digraph->aGCD) {
					int ibf0 = digraph->ibf2(t);
					int ibf1 = digraph->ibf2_fast(t);

					if (ibf0 != ibf1) {
						cerr << "At time " << t << " the two algorithms have different values: " << ibf0 << "\t" << ibf1 << endl;
						exit(EXIT_FAILURE);
					}
				}
			}
			*/
		}
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis(digraphs,i,linearization);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphIBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis(Digraph** digraphs, int i, bool linearization) {
	Digraph* digraph = digraphs[i];
	for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
		Node* node = *iter;
		bool schedulable = digraph_ibf_analysis(digraphs, i, linearization, node->wcet, node->deadline);
		if (!schedulable) {
			cout << "Digraph-" << i << ": (wcet,deadline)=<" <<node->wcet<< "," << node->deadline << ">" << endl; 	
			return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis(Digraph** digraphs, int i, bool linearization, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->aGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			if (linearization)
				sum += digraphj->ibf(t);
			else
				sum += digraphj->ibf2_fast(t);
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return false;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice,bool periodicityProperty, bool deadlineProperty) {
	Timer timer;
	timer.start();

	if (deadlineProperty) { // the l-Mad property
		prepare_all_digraphs(digraphs,n,2);

		for (int i=0; i<n; i++) {
			Digraph* digraph = digraphs[i];

			//cout << "digraph-"<< i << "=>tf=" << digraph->tf2 << endl;

			if (periodicityProperty && digraph->strongly_connected) {
				digraph->calculate_linear_factor();
				digraph->unit_digraph = new UnitDigraph(digraph);
				digraph->unit_digraph->prepare3(false);
			}
		}
	}
	else {
		for (int i=0; i<n; i++) {
			Digraph* digraph = digraphs[i];
			digraph->calculate_period_gcd();
			digraph->calculate_all_gcd();

			if (periodicityProperty) {
				digraph->generate_strongly_connected_components();
				digraph->check_strongly_connected();

				if (digraph->strongly_connected) {
					digraph->calculate_linear_factor();
					digraph->unit_digraph = new UnitDigraph(digraph);
					digraph->unit_digraph->prepare3(false);
				}
			}
		}
	}

	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		nDigraphNode += digraph->node_vec.size();
		nDigraphEdge += digraph->edge_vec.size();

		if (deadlineProperty) {
			tf_max = max(tf_max,digraph->tf2);
		}
		else {
			for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
				tf_max = max(tf_max,(*iter)->deadline);
			}
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;
		int tf = ceil(1.0*tf_max/digraph->aGCD)*digraph->aGCD;

		if (periodicityProperty && digraph->strongly_connected) {
			nIrreducibleDigraphs ++;
			digraph->calculate_ibf_with_periodicity_DP(tf,digraph->ibf_map2_DP);
		}
		else {
			digraph->calculate_ibf_without_periodicity_DP(tf,digraph->ibf_map2_DP);	
		}

		if (false) {
			cout << "==============" << "Digraph-" << i << "==============" << endl;
			for (map<int,int>::iterator mIter = digraph->ibf_map2_DP.begin(); mIter != digraph->ibf_map2_DP.end(); mIter++) 
				cout << "Digraph-" << i << ".ibf(" << mIter->first << ")=" <<mIter->second << endl;
		}

		if (deadlineProperty) digraph->calculate_dbf_without_periodicity_DP(digraph->tf2,digraph->dbf_map2_DP);
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(digraphs,i,deadlineProperty);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphIBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i,bool deadlineProperty) {
	Digraph* digraph = digraphs[i];
	if (deadlineProperty) {
		int gcd = 0;
		for (int j=0; j<=i; j++) {
			Digraph* digraphj = digraphs[j];
			gcd = Utility::math_gcd(gcd, digraphj->aGCD);
		}

		int tf = digraph->tf2/digraph->pGCD;

		for (int prime = 0;  prime<= tf; prime++) {
			int t = prime * digraph->pGCD;

			bool schedulable = false;
			for (int tdot = 0; tdot <= t; tdot += gcd) {
				double sum = digraph->dbf2(t,digraph->dbf_map2_DP);
				for (int j=0; j<i; j++) {
					Digraph* digraphj = digraphs[j];
			
					sum += digraphj->ibf2(tdot,digraphj->ibf_map2_DP);
				}
				if (sum <= tdot) {
					schedulable = true;
					break;
				}
				tdot = floor(1.0*sum/gcd)*gcd;
			}

			if (!schedulable) {
				cout << "Digraph-" << i << ": prime = " <<prime << endl;
				return false;
			}
		}
		return true;
	}
	else {
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			Node* node = *iter;
			bool schedulable = digraph_ibf_analysis_dynamic_programming(digraphs, i, node->wcet, node->deadline);
			if (!schedulable) {
				cout << "Digraph-" << i << ": (wcet,deadline)=<" <<node->wcet<< "," << node->deadline << ">" << endl; 	
				return false;
			}
		}
		return true;
	}
}

bool SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->aGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->ibf2(t,digraphj->ibf_map2_DP);
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return false;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int n, int choice,bool periodicityProperty) {
	Timer timer;
	timer.start();

	prepare_all_digraphs(digraphs,n,2);

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		if (periodicityProperty && digraph->strongly_connected) {
			digraph->unit_digraph = new UnitDigraph(digraph);
			digraph->unit_digraph->prepare3(false);
		}
	}

	// calculate \mathcal{L} and the max response time
	double maxL = 0;
	//double maxRT = 0;
	int maxDeadline = INT_MIN;

	for (vector<Node*>::iterator iter = digraphs[n-1]->node_vec.begin(); iter != digraphs[n-1]->node_vec.end(); iter++) {
		maxDeadline = max(maxDeadline, (*iter)->deadline);
	}

	double sumU = 0;
	double sumC1 = 0;
	//double sumC2 = 0;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		sumU += digraph->linear_factor;
		sumC1 += digraph->c_ibf;
		//if (i!=n-1) sumC2 += digraph->c_ibf;
	}
	maxL = sumC1 / (1-sumU);
	//maxRT = (sumC2 + node->wcet) / (1-sumU);

	//cout << "maxL = " << maxL << ", maxDeadline = " << maxDeadline << ", linear factor = " << digraphs[n-1]->linear_factor << endl;

	for (int i=0; i<n-1; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = maxL+maxDeadline;
		int tf = ceil((maxL+maxDeadline)/digraph->aGCD)*digraph->aGCD;
		//cout << "maximum length is " << tf << endl;
		//cout << "Preparing digraph-" << i << endl;
		int maxSeperation = INT_MIN;
		if (periodicityProperty) {
			for (auto edge : digraph->edge_vec)
				maxSeperation = max(maxSeperation,edge->separationTime);
		}

		if (periodicityProperty && digraph->strongly_connected /* && 1.0*tf/maxSeperation >= 50 */) {
			nIrreducibleDigraphs ++;
			digraph->calculate_ibf_with_periodicity_DP(tf,digraph->ibf_map2_DP);
		}
		else {
			digraph->calculate_ibf_without_periodicity_DP(tf,digraph->ibf_map2_DP);	
		}

		if (false) {
			cout << "==============" << "Digraph-" << i << "==============" << endl;
			for (map<int,int>::iterator mIter = digraph->ibf_map2_DP.begin(); mIter != digraph->ibf_map2_DP.end(); mIter++) 
				cout << "Digraph-" << i << ".ibf(" << mIter->first << ")=" <<mIter->second << endl;
		}

	}

	for (int i=n-1; i<n; i++) {
		bool schedulable = SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs,i);
		if (!schedulable) {
			timer.end();
			tDigraphIBF += timer.getTime();
			return false;
		}
	}
	
	//schedulable = SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs,n-1,maxL, maxRT,node);

	timer.end();
	tDigraphIBF += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i) {
	Digraph* digraphi = digraphs[i];

	double maxL = 0;

	double sumU = 0;
	double sumC = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];

		sumU += digraphj->linear_factor;
		sumC += digraphj->c_ibf;
	}
	maxL = sumC / (1-sumU);

	for (vector<Node*>::iterator iter = digraphi->node_vec.begin(); iter != digraphi->node_vec.end(); iter++) {
		Node* node = *iter;
		bool schedulable = digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs, i, maxL, node);
		if (!schedulable) {
			//cout << "Digraph-" << i << ": (wcet,deadline)=<" <<node->wcet<< "," << node->deadline << ">" << endl; 	
			return false;
		}
	}

	return true;
}

bool SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i, double maxL, Node* node) {
	// initial rts
	set<RequestTriple*> rts; // used to record the valid paths ending at the node
	RequestTriple* rt0 = new RequestTriple(node->wcet, 0, node);
	rts.insert(rt0);

	set<RequestTriple*> all_rts;
	all_rts.insert(rt0);

	map<int,int> mPair; // used to record the request triples which have been tested

	// check the path containing only itself
	bool schedulable = digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs,i,maxL,node->wcet,node->deadline,0);
	if (!schedulable) {
		//cout << "++++++++++++++" << endl;
		goto END;
	}

	// check the path containing other vertices
	while (!rts.empty()) {
		set<RequestTriple*> temp_rts;

		for (set<RequestTriple*>::iterator iter = rts.begin(); iter != rts.end(); iter++) {
			RequestTriple* rt = *iter;
			Node* snk = rt->node;

			for (list<Edge*>::iterator eIter = snk->in.begin(); eIter != snk->in.end(); eIter++) {
				Edge* edge = *eIter;
				Node* src = edge->src_node;
				int e = rt->e + src->wcet;
				int r = rt->r + edge->separationTime;

				if (r > maxL) continue; 

				bool hasNext = true;

				// check whether the prefix path (without the last vertex) can interfere the execution of the last vertex
				double sum = e-node->wcet;
				for (int j=0; j<i; j++) {
					Digraph* digraphj = digraphs[j];

					sum += digraphj->ibf2(r,digraphj->ibf_map2_DP);
				}
				if (sum <= r) hasNext = false;

				if (mPair.find(r)!=mPair.end()) {
					if (mPair[r] >= e)
						schedulable = true;
				}
				else {
					schedulable = digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs,i,maxL,e,node->deadline,r);
					if (!schedulable) goto END;
					mPair[r] = e;
				}

				if (hasNext) {
					RequestTriple* nrt = new RequestTriple(e,r,src);
					temp_rts.insert(nrt);
					all_rts.insert(nrt);
				}
			}
		}
		rts = temp_rts;
	}

END:
	// release rts
	for (set<RequestTriple*>::iterator iter = all_rts.begin(); iter != all_rts.end(); iter++)
		delete *iter;

	return schedulable;
}

double SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i, double maxL, RequestTriple* rt) {
	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->aGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=2*maxL; t+=gcd) {
		double sum = rt->e;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->ibf2(t,digraphj->ibf_map2_DP);
		}
		if (sum <= t) return sum-rt->r;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return 0;
}

bool SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(Digraph** digraphs, int i, double maxL, int wcet,int deadline, int releaseTime) {
	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->aGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=releaseTime+deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->ibf2(t,digraphj->ibf_map2_DP);
		}
		if (sum <= t) {
			//cout << "<" << wcet << "," << deadline << "," << releaseTime << ", =>>> rt = " << sum - releaseTime << endl;
			return true;
		}

		t = floor(1.0*sum/gcd)*gcd;
	}
	return false;
}

bool SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_constrained_deadline(Digraph** digraphs, int n, int choice,bool periodicityProperty) {
	Timer timer;
	timer.start();

	prepare_all_digraphs(digraphs,n,2);

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];

		if (periodicityProperty && digraph->strongly_connected) {
			digraph->unit_digraph = new UnitDigraph(digraph);
			digraph->unit_digraph->prepare3(false);
		}
	}

	Digraph* lowDigraph = digraphs[n-1];
	Node* node = lowDigraph->node_vec.at(0);
	for (vector<Node*>::iterator iter = lowDigraph->node_vec.begin(); iter != lowDigraph->node_vec.end(); iter++) {
		if ((*iter)->wcet > node->wcet) node = *iter;
	}

	double sumU = 0;
	double sumC = 0;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		sumU += digraph->linear_factor;
		if (i!=n-1) sumC += digraph->c_ibf;
	}
	// calculate the max response time
	double maxRT = (sumC + node->wcet) / (1-sumU);

	//cout << "maxL = " << maxL << ", maxRT = " << maxRT << endl;

	for (int i=0; i<n-1; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = maxRT;
		int tf = ceil((maxRT)/digraph->aGCD)*digraph->aGCD;

		if (periodicityProperty && digraph->strongly_connected) {
			nIrreducibleDigraphs ++;
			digraph->calculate_ibf_with_periodicity_DP(tf,digraph->ibf_map2_DP);
		}
		else {
			digraph->calculate_ibf_without_periodicity_DP(tf,digraph->ibf_map2_DP);	
		}

		if (false) {
			cout << "==============" << "Digraph-" << i << "==============" << endl;
			for (map<int,int>::iterator mIter = digraph->ibf_map2_DP.begin(); mIter != digraph->ibf_map2_DP.end(); mIter++) 
				cout << "Digraph-" << i << ".ibf(" << mIter->first << ")=" <<mIter->second << endl;
		}

	}

	bool schedulable = false;
	/*
	for (int i=n-1; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_arbitrary_deadline(digraphs,i,maxL);
		if (!schedulable) break;
	}
	*/

	double nodeRT = SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_constrained_deadline(digraphs,n-1,node->wcet, maxRT);

	//cout << "wcet = " << node->wcet << ", ibf_rt = " << nodeRT << ", lubibf_RT = " << maxRT << endl; 

	timer.end();
	tDigraphIBF += timer.getTime();
	return schedulable;
}

double SchedulabilityAnalysis::digraph_rt_ibf_analysis_dynamic_programming_under_constrained_deadline(Digraph** digraphs, int i, int wcet, double maxRT) {
	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->aGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=maxRT; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->ibf2(t,digraphj->ibf_map2_DP);
		}
		if (sum <= t) return sum;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return 0;
}

bool SchedulabilityAnalysis::reachability_digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
	}
	
	generate_critical_vertices_for_reachability_digraphs(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		//nDigraphNode += digraph->node_vec.size();
		//nDigraphEdge += digraph->edge_vec.size();

		for (vector<Edge*>::iterator iter = digraph->cedge_vec.begin(); iter != digraph->cedge_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;
		
		digraph->calculate_ibf_without_periodicity_on_reachability_digraph_DP(tf_max,digraph->ibf_map2_DP);	

		if (false) {
			cout << "==============" << "Digraph-" << i << "==============" << endl;
			for (map<int,int>::iterator mIter = digraph->ibf_map2_DP.begin(); mIter != digraph->ibf_map2_DP.end(); mIter++) 
				cout << "Digraph-" << i << ".ibf(" << mIter->first << ")=" <<mIter->second << endl;
		}
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::reachability_digraph_ibf_analysis_dynamic_programming(digraphs,i);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphIBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::reachability_digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i) {
	Digraph* digraph = digraphs[i];
	for (vector<Edge*>::iterator iter = digraph->cedge_vec.begin(); iter != digraph->cedge_vec.end(); iter++) {
		Edge* edge = *iter;
		bool schedulable = reachability_digraph_ibf_analysis_dynamic_programming(digraphs, i, edge->weight, edge->deadline);
		if (!schedulable) {
			cout << "Digraph-" << i << ": (wcet,deadline)=<" <<edge->weight<< "," << edge->deadline << ">" << endl; 	
			return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::reachability_digraph_ibf_analysis_dynamic_programming(Digraph** digraphs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->aGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->ibf2(t,digraphj->ibf_map2_DP);
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return false;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis_with_rbf_linear_defect(Digraph** digraphs, int n, int choice) {
	Timer timer;
	timer.start();

	prepare_all_digraphs(digraphs,n,choice);
	
	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->tf = tf_max;

		digraph->calculate_ibf_with_periodicity_by_rbf_linear_defect(tf_max, digraph->ibf_map2_fast);
	}

	bool schedulable = false;
	for (int i=0; i<n; i++) {
		schedulable = SchedulabilityAnalysis::digraph_ibf_analysis_with_rbf_linear_defect(digraphs,i);
		if (!schedulable) break;
	}

	timer.end();
	tDigraphIBF += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis_with_rbf_linear_defect(Digraph** digraphs, int i) {
	Digraph* digraph = digraphs[i];
	for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
		Node* node = *iter;
		bool schedulable = digraph_ibf_analysis_with_rbf_linear_defect(digraphs, i, node->wcet, node->deadline);
		if (!schedulable) {
			cout << "Digraph-" << i << ": (wcet,deadline)=<" <<node->wcet<< "," << node->deadline << ">" << endl; 	
			return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::digraph_ibf_analysis_with_rbf_linear_defect(Digraph** digraphs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Digraph* digraphj = digraphs[j];
		gcd = Utility::math_gcd(gcd,digraphj->aGCD);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Digraph* digraphj = digraphs[j];
			
			sum += digraphj->ibf2_fast(t);
			
		}
		if (sum <= t) return true;

		t = floor(1.0*sum/gcd)*gcd;
	}
	return false;
}

void SchedulabilityAnalysis::prepare_all_stateflows(Stateflow** sfs, int n, int choice) {
	Timer timer;
	// prepare linear upper bounds and the maximal time length
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		//sf->write_graphviz(cout);
		//sf->precise_digraph->write_graphviz(cout);

		//sf->check_irreducible();

		if (choice == 0) {
			sf->generate_rbfs();
			sf->calculate_linear_factor();
			timer.start();
			sf->calculate_csum();
			timer.end();
			tCalCSum += timer.getTime();

			timer.start();
			sf->calculate_tf0(sfs,i);
			timer.end();
			tCalTF0 += timer.getTime();
		}
		// cout << "isIrred=" << sf->isIrred <<endl;
		// sf->calculate_generialized_period();
		// sf->calculate_generialized_defect();
		
		//sf->generate_simple_digraph();

		if (choice == 1 || choice == 2) {
			timer.start();
			sf->generate_precise_digraph();
			sf->calculate_linear_upper_bounds(true);
			timer.end();
			tCalLinearBounds += timer.getTime();

			timer.start();
			sf->calculate_tf1(sfs,i);
			timer.end();
			tCalTF1 += timer.getTime();

			timer.start();
			sf->calculate_tf2(sfs,i);
			timer.end();
			tCalTF2 += timer.getTime();

			bpRatio0 += sf->tf1/sf->tf0;
			bpRatio1 += sf->tf2/sf->tf0;
		}
		//cout<<"tf0="<<sf->tf0<<"\ttf1="<<sf->tf1<<"\ttf2="<<sf->tf2<<endl;
	}

	// generate all execution request matrices
	//calculate_request_execution_matrices_without_periodicity_property(sfs,n,choice);
	cout << "Generating critical action pairs" <<endl;
	generate_critical_action_pair(sfs,n);
	if (false) output_critical_action_pair(sfs,n,cout);
	cout << "size of ActionPair = " << sizeof(ActionPair) << endl;
}

void SchedulabilityAnalysis::prepare_all_stateflows2(Stateflow** sfs, int n, bool RBF, bool periodicity) {
	Timer timer;
	timer.start();
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		//sf->write_graphviz(cout);
		sf->calculate_gcd();
		sf->calculate_t_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		if (RBF)
			sf->generate_rbfs();
		else
			sf->generate_ibfs();
		// show rbfs
		if (false) {
			for (int i=0; i<sf->n_state; i++) for (int j=0; j<sf->n_state; j++) {
				cout<<"rbf_{"<<i<<","<<j<<"}="<<endl;
				Utility::output_matrix(sf->rbfs[i][j], sf->n_time_instance,sf->n_time_instance);
			}
		}

		// statistics the average degree for stateflows
		avgDegree += 1.0*sf->trans.size()/sf->states.size();

		// show execution request matrix
		if (false) Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);

		if (periodicity) {
			sf->calculate_linear_factor();
		    //cout << sf->isIrred << endl;
			if (sf->isIrred) {
				sf->calculate_generialized_period();
			}
		}
		else 
			sf->isIrred = false; // close the calculation of the periodicity parameters
		
	}

	//cout << "Generating critical action pairs" <<endl;
	generate_critical_action_pair(sfs,n);
	if (false) output_critical_action_pair(sfs,n,cout);
	//cout << "size of ActionPair = " << sizeof(ActionPair) << endl;

	int tf_max = INT_MIN;

	int sum = 0;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		for (map<int,vector<ActionPair>>::iterator mIter = sf->mCAP.begin(); mIter != sf->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;

			for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
				sum += sizeof(ActionPair);
				ActionPair ap = *apIter;
				tf_max = max(tf_max,ap.stime+ap.d);
			}
		}
	}
	//cout << "sum size of ActionPair = " << sum << "\tt_f=" << tf_max << endl;

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->tf0 = tf_max/sf->hyperperiod+2;
		if (periodicity && sf->isIrred) {
			sf->calculate_generialized_defect();
			//if (sf->gdef != -1) cout << "ldef=" << sf->gdef << "\ttf=" << sf->tf0 << endl;
		}

		if (false) cout <<"matrix power=" << tf_max/sf->hyperperiod+2 << endl;

		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+2);
	}

	//prepare_all_stateflows(sfs,n,choice);

	timer.end();
	tPrepareStateflow += timer.getTime();
}

void SchedulabilityAnalysis::reset_calculating_containers(Stateflow** sfs, int n) {
	// reset all the calculating containers
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		
		// release rbf_vec
		for (vector<StartFinishTime*>::iterator iter = sf->rbf_vec.begin(); iter != sf->rbf_vec.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
		sf->rbf_vec.clear();

		// release ibf_vec
		for (vector<StartFinishTime*>::iterator iter = sf->ibf_vec.begin(); iter != sf->ibf_vec.end(); iter++) {
			delete *iter;
			*iter = NULL;
		}
		sf->ibf_vec.clear();
	}
}

void SchedulabilityAnalysis::calculate_request_execution_matrices_without_periodicity_property(Stateflow** sfs, int n, int choice) {
	Timer timer;
	//timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}
	timer.start();
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		if (true) cout <<"matrix power=" << tf_max/sf->hyperperiod+2 << endl;
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+2);
	}
	timer.end();
	tCalRequestExecutionMatrixWithoutPeriodicityProperty += timer.getTime();
}

void SchedulabilityAnalysis::calculate_request_execution_matrices_with_periodicity_property(Stateflow** sfs, int n, int choice) {
	Timer timer;
	//timer.start();
	// set tf, denoted as the maximum time length for the special stateflow
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		int tf;
		if (choice == 0) { 
			if (sf->tf0 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf0);
		}
		else if (choice == 1) {
			if (sf->tf1 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf1);
		}
		else if (choice == 2) {
			if (sf->tf2 == POS_INFINITY) {
				cerr << "Should not arrive here!" <<endl;
				exit(EXIT_FAILURE);
			}
			tf = ceil(sf->tf2);
		}
		else { 
			cerr << "Error choice=" << choice <<endl;
			exit(EXIT_FAILURE);
		}
		tf_max = max(tf_max,tf);
	}
	timer.start();
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		sf->calculate_exec_req_matrix_power(tf_max/sf->hyperperiod+2);
	}
	timer.end();
	tCalRequestExecutionMatrixWithPeriodicityProperty += timer.getTime();
}

void SchedulabilityAnalysis::generate_critical_action_pair(Stateflow** sfs, int n) {
	Timer timer;
	timer.start();
	int hyperperiod = 1;
	typedef vector<Transition*>::iterator TranVecIter;
	typedef list<Transition*>::iterator TranListIter;
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		hyperperiod = Utility::math_lcm(hyperperiod, sf->hyperperiod);
		int maxDeadline = INT_MIN;
		// generate all the action pairs
		for(TranVecIter iter = sf->trans.begin(); iter != sf->trans.end(); iter ++) {
			Transition* tran = *iter;
			// generate all the action pairs for this transition
			for (int stime = 0; stime < hyperperiod; stime += tran->period) {
				State* snk = tran->snk;
				int deadline = INT_MAX;

				if (snk->out.empty()) deadline = stime+tran->period;

				for (TranListIter iter2 = snk->out.begin(); iter2 != snk->out.end(); iter2++) {
					int t = 0;
					while (t<=stime) t+=(*iter2)->period; 
					deadline = min(deadline,t);
				}
				deadline = deadline - stime;
				maxDeadline = max(maxDeadline, deadline);
				ActionPair ap = ActionPair(tran->wcet,deadline,stime,i,tran);
				vector<ActionPair>& vec_AP = sf->mCAP[stime];
				// Determine intra-FSM dominances
				int nSize = vec_AP.size();
				bool* found = new bool[nSize];
				for (int k=0; k<nSize; k++) found[k] = false;
				bool dominated = false;
				bool flag = false;
				int nItem = 0;
				for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
					ActionPair nap = *apIter;
					if (nap>ap) { dominated = true; break;}
					if (ap>nap) { flag = true; found[nItem] = true;}
					nItem ++;
				}

				if (flag) {
					vector<ActionPair> tempAP;
					nItem = 0;
					for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
						ActionPair nap = *apIter;
						if (!found[nItem]) tempAP.push_back(nap); 
						nItem ++;
					}
					tempAP.push_back(ap);
					sf->mCAP[stime] = tempAP;
				}
				else {
					if (!dominated)	
						vec_AP.push_back(ap);
				}

				delete[] found;
			}
		}

		sf->maxDeadline = maxDeadline;

		// Determine inter-FSM dominances
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			map<int,vector<ActionPair>>::iterator mIter;
			for(mIter = sfj->mCAP.begin(); mIter != sfj->mCAP.end(); mIter++) {
				int key = (*mIter).first;
				vector<ActionPair> vec_AP = (*mIter).second;

				if (sf->mCAP.find(key) == sf->mCAP.end()) continue;

				vector<ActionPair> vec_AP2 = sf->mCAP[key];
				int nSize = vec_AP.size();
				bool* found = new bool[nSize];
				for (int k=0; k<nSize; k++) found[k] = false;
				bool flag = false;
				int nItem = 0;

				for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
					ActionPair ap = *apIter;

					for (vector<ActionPair>::iterator apIter2 = vec_AP2.begin(); apIter2 != vec_AP2.end(); apIter2++) {
						ActionPair ap2 = *apIter2;
						if (ap2>ap) { flag = true; found[nItem] = true; break;}
					}
					nItem ++;
				}

				if (flag) {
					vector<ActionPair> tempAP;
					nItem = 0;

					for (vector<ActionPair>::iterator apIter = vec_AP.begin(); apIter != vec_AP.end(); apIter++) {
						ActionPair ap = *apIter;

						if (!found[nItem]) tempAP.push_back(ap);
						nItem++;
					}

					sfj->mCAP[key] = tempAP;
				}
				delete[] found;
			}
		}
	}
	timer.end();
	tGenerateCriticalActionPairs += timer.getTime();
}

bool SchedulabilityAnalysis::generate_request_function(Stateflow** sfs, int n, bool output) {
	Timer timer;
	timer.start();
	for (int i=0; i<n-1; i++) {
		Stateflow* sfi = sfs[i];
		int maxDeadline = INT_MIN;
		int tk = INT_MIN;
		for (int j=i+1; j<n; j++) {
			Stateflow* sfj = sfs[j];
			maxDeadline = max(maxDeadline, sfj->maxDeadline);

			for (map<int,vector<ActionPair>>::iterator mIter = sfj->mCAP.begin(); mIter != sfj->mCAP.end(); mIter++) {
				int first = mIter->first;
				tk = max(tk,first);
			}
		}
		// generate request functions
		bool success = generate_request_function(sfi,i,tk+maxDeadline,output);
		if (!success) {
			timer.end();
			tGenerateRequestFunctions+=timer.getTime();
			return false;
		}
	}

	timer.end();
	tGenerateRequestFunctions+=timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::generate_request_function(Stateflow* sf, int i, int maxDeadline, bool output) {
	if (maxDeadline == 0) maxDeadline = sf->hyperperiod;
	int gcd = sf->gcd;
	int hyperperiod = sf->hyperperiod;
	vector<RequestFunction> vec_rf;

	int t = 0;
	int count = 0;

	while (t <= hyperperiod + maxDeadline) {
		for (set<int>::iterator iter = sf->rbf_time_instances.begin(); iter != sf->rbf_time_instances.end(); iter++) {
			int instance = *iter;
			if (instance == hyperperiod) continue;

			t = count*hyperperiod + instance;
			if (t > hyperperiod+maxDeadline) break;

			vector<RequestFunction> temp_vec_rf;
			if (t==0) {
				for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
					Transition* tran = *iter;
					if (t%tran->period!=0) continue;
					RequestFunction rf = RequestFunction(i,gcd,hyperperiod,t,t,hyperperiod+maxDeadline);
					rf.update_next_action(t,tran);
					temp_vec_rf.push_back(rf);
				}
				vec_rf.clear();
				vec_rf = temp_vec_rf;
				temp_vec_rf.clear();
				if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
				continue;
			}

			for (vector<RequestFunction>::iterator rfIter = vec_rf.begin(); rfIter != vec_rf.end(); rfIter++) {
				RequestFunction rf = *rfIter;
				list<Transition*> list_tran = rf.lastTran->snk->out;
				list<Transition*>::iterator tranIter;
				for (tranIter = list_tran.begin(); tranIter != list_tran.end(); tranIter++) {
					RequestFunction temp_rf2 = RequestFunction(rf);
					Transition* tran = *tranIter;
					if (t%tran->period==0) {
						temp_rf2.update_next_action(t,tran);
					} else {
						temp_rf2.update_next_action(t,NULL);
					}
					temp_vec_rf.push_back(temp_rf2);
				}

				if (list_tran.empty()) {
					rf.update_next_action(t, NULL);
					temp_vec_rf.push_back(rf);
				}
			}

			// add other request functions triggered within one hyperperiod but 0
			if (t < hyperperiod) {
				for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
					Transition* tran = *iter;
					if (t%tran->period!=0) continue;
					RequestFunction rf = RequestFunction(i,gcd,hyperperiod,t,t,hyperperiod+maxDeadline);
					rf.update_next_action(t,tran);
					temp_vec_rf.push_back(rf);
				}
			}

			if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
			if (false) {
				output_critical_request_function(temp_vec_rf,cout);
				cout << "=========================" <<endl;
			}
			intra_loop_dominance(temp_vec_rf);
			if (false) {
				output_critical_request_function(temp_vec_rf,cout);
				cout << "+++++++++++++++++++++++++" << endl;
			}
			vec_rf.clear();
			vec_rf = temp_vec_rf;
			temp_vec_rf.clear();
			if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
			//if (vec_rf.size() >= 500) return false;
		}
		count++;
	}

	//output_critical_request_function(vec_rf,cout);
	end_loop_dominance(vec_rf);
	//output_critical_request_function(vec_rf,cout);
	sf->CRF = vec_rf;
	cout << "upper bound=" << maxDeadline  << "\tnumber of request functions = " << vec_rf.size() << endl;
	return true;
}

void SchedulabilityAnalysis::generate_request_function_for_each_action_instance(Stateflow** sfs, int n, bool output) {
	Timer timer;
	timer.start();
	for (int i=1; i<n; i++) {
		Stateflow* sfi = sfs[i];

		for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
			int maxDeadline = INT_MIN;
			vector<ActionPair> aps = mIter->second;
			if (aps.empty()) continue;

			for (vector<ActionPair>::iterator iter = aps.begin(); iter != aps.end(); iter++) {
				maxDeadline = max(maxDeadline,iter->d);
			}

			int tk = mIter->first;
			int n = tk/sfi->hyperperiod;
			int r = tk%sfi->hyperperiod;
			if (n==0 && r==0) tk=0;
			else if (n != 0 && r ==0) {
				tk=sfi->hyperperiod;
				int temp = *(--sfi->rbf_time_instances.lower_bound(tk));
				tk = tk-temp;
			} 
			else {
				tk = r;
				int temp = *(--sfi->rbf_time_instances.lower_bound(tk));
				tk = tk-temp;
			} 

			if (mIter->first-tk < 0) {
				cerr << "Error start time" << mIter->first-tk <<endl;
				exit(EXIT_FAILURE);
			}
			
			generate_request_function_for_each_action_instance(sfs, i, mIter->first-tk,mIter->first+maxDeadline,output);
		}
	}
	timer.end();
	tGenerateRequestFunctions+=timer.getTime();
}

void SchedulabilityAnalysis::generate_request_function_for_each_action_instance(Stateflow** sfs, int i, int tk, int maxDeadline, bool output) {
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		generate_request_function_for_each_action_instance(sfj,j,tk,maxDeadline,output);
	}
}

void SchedulabilityAnalysis::generate_request_function_for_each_action_instance(Stateflow* sf, int i, int tk, int maxDeadline, bool output) {
	if (maxDeadline == 0) maxDeadline = sf->hyperperiod;
	int gcd = sf->gcd;
	int hyperperiod = sf->hyperperiod;

	vector<RequestFunction> vec_rf;

	int t = 0;
	int count = 0;
	bool flag = true;

	while (t <= maxDeadline) {
		for (set<int>::iterator iter = sf->rbf_time_instances.begin(); iter != sf->rbf_time_instances.end(); iter++) {
			int instance = *iter;
			if (instance == hyperperiod) continue;

			t = count*hyperperiod + instance;
			if (t < tk) continue;
			if (t > maxDeadline) break;

			vector<RequestFunction> temp_vec_rf;

			if (flag) {
				for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
					Transition* tran = *iter;
					if (t%tran->period!=0) continue;
					RequestFunction rf = RequestFunction(i,gcd,hyperperiod,t,t,hyperperiod+maxDeadline);
					rf.update_next_action(t,tran);
					temp_vec_rf.push_back(rf);
				}
				vec_rf.clear();
				vec_rf = temp_vec_rf;
				temp_vec_rf.clear();
				if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
				flag = false;
				continue;
			}

			for (vector<RequestFunction>::iterator rfIter = vec_rf.begin(); rfIter != vec_rf.end(); rfIter++) {
				RequestFunction rf = *rfIter;
				list<Transition*> list_tran = rf.lastTran->snk->out;
				list<Transition*>::iterator tranIter;
				for (tranIter = list_tran.begin(); tranIter != list_tran.end(); tranIter++) {
					RequestFunction temp_rf2 = RequestFunction(rf);
					Transition* tran = *tranIter;
					if (t%tran->period==0) {
						temp_rf2.update_next_action(t,tran);
					} else {
						temp_rf2.update_next_action(t,NULL);
					}
					temp_vec_rf.push_back(temp_rf2);
				}

				if (list_tran.empty()) {
					rf.update_next_action(t, NULL);
					temp_vec_rf.push_back(rf);
				}
			}

			// add other request functions triggered within one hyperperiod but 0
			if (t < hyperperiod) {
				for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
					Transition* tran = *iter;
					if (t%tran->period!=0) continue;
					RequestFunction rf = RequestFunction(i,gcd,hyperperiod,t,t,hyperperiod+maxDeadline);
					rf.update_next_action(t,tran);
					temp_vec_rf.push_back(rf);
				}
			}

			if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
			if (false) {
				output_critical_request_function(temp_vec_rf,cout);
				cout << "=========================" <<endl;
			}
			intra_loop_dominance(temp_vec_rf);
			if (false) {
				output_critical_request_function(temp_vec_rf,cout);
				cout << "+++++++++++++++++++++++++" << endl;
			}
			vec_rf.clear();
			vec_rf = temp_vec_rf;
			temp_vec_rf.clear();
			if (output) cout << i << "=>>>>>>" << vec_rf.size() <<endl;
			//if (vec_rf.size() >= 500) return false;
		}
		count++;
	}

	for (vector<RequestFunction>::iterator iter = vec_rf.begin(); iter != vec_rf.end(); iter++) {
		RequestFunction& rf = *iter;
		if (rf.finish < maxDeadline)
			rf.update_next_action(maxDeadline,NULL);
	}

	/*
	if (vec_rf.empty()) {
		cerr << "There is no request function!" << endl;
		exit(EXIT_FAILURE);
	}
	*/

	//output_critical_request_function(vec_rf,cout);
	if (vec_rf.empty()) {
		RequestFunction rf = RequestFunction(i,gcd,hyperperiod,tk,tk,hyperperiod+maxDeadline);
		vec_rf.push_back(rf);
	}
	
	end_loop_dominance(vec_rf);

	//output_critical_request_function(vec_rf,cout);
	sf->mmCRF[tk][maxDeadline] = vec_rf;
	if (output) 
		cout << "Stateflow" << i << ".[" << tk << "," << maxDeadline << "]"  << "\tnumber of request functions = " << vec_rf.size() << endl;
}

void SchedulabilityAnalysis::intra_loop_dominance(vector<RequestFunction>& vec_rf) {
	// detemine the dominance for rfs
	int nSize = vec_rf.size();
	vector<bool> found = vector<bool>(nSize,false);
	vector<int> count = vector<int>(nSize,0);

	for (int j=0; j<nSize; j++) {
		RequestFunction rf1 = vec_rf.at(j);
		if (found[j]) continue;
		for (int k=0; k<nSize; k++) {
			if (k==j) continue;
			if (found[k]) continue;
			RequestFunction rf2 = vec_rf.at(k);
			if (rf1.start == rf2.start && rf1.lastTran == rf2.lastTran && rf1.isCounted == rf2.isCounted && rf1 > rf2) { 
				found[k] = true;
				count[j]++;
			}
		}
	}

	vector<RequestFunction> temp_vec_rf;
	for (int k=0; k<nSize; k++) {
		if (!found[k]) temp_vec_rf.push_back(vec_rf.at(k));
	}

	if (temp_vec_rf.empty()) {
		int index = 0;
		int temp = INT_MIN;
		for (int k=0; k<nSize; k++) {
			if (count[k] > temp) {
				temp = count[k];
				index = k;
			}
		}
		//cout << vec_rf.size() << "\t" << index <<endl;
		temp_vec_rf.push_back(vec_rf.at(index));
	}

	vec_rf.clear();
	vec_rf = temp_vec_rf;
	temp_vec_rf.clear();

	found.clear();
	count.clear();
}

void SchedulabilityAnalysis::end_loop_dominance(vector<RequestFunction>& vec_rf) {
	// detemine the dominance for rfs
	int nSize = vec_rf.size();
	vector<bool> found = vector<bool>(nSize,false);
	vector<int> count = vector<int>(nSize,0);
	

	for (int j=0; j<nSize; j++) {
		if (found[j]) continue;
		RequestFunction rf1 = vec_rf.at(j);
		for (int k=0; k<nSize; k++) {
			if (k==j) continue;
			if (found[k]) continue;
			RequestFunction rf2 = vec_rf.at(k);
			if (rf1.start == rf2.start && rf1 > rf2) { 
				found[k] = true;
				count[j]++;
			}
		}
	}

	vector<RequestFunction> temp_vec_rf;
	for (int k=0; k<nSize; k++) {
		if (!found[k]) temp_vec_rf.push_back(vec_rf.at(k));
	}

	if (temp_vec_rf.empty()) {
		int index = 0;
		int temp = INT_MIN;
		for (int k=0; k<nSize; k++) {
			if (count[k] > temp) {
				temp = count[k];
				index = k;
			}
		}
		temp_vec_rf.push_back(vec_rf.at(index));
	}

	vec_rf.clear();
	vec_rf = temp_vec_rf;
	temp_vec_rf.clear();

	found.clear();
	count.clear();
}


void SchedulabilityAnalysis::generate_request_function_abstract_tree(Stateflow** sfs, int n, bool output) {
	Timer timer;
	timer.start();
	for (int i=1; i<n; i++) {
		Stateflow* sfi = sfs[i];
		for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;
			for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
				ActionPair ap = *apIter;
				int tk = 0;
				for (int align = 0; align <=ap.stime; align += sfi->hyperperiod) {
					bool check = true;
					for (int j=0; j<i; j++) {
						Stateflow* sfj = sfs[j];
						if (align % sfj->hyperperiod != 0) {
							check = false;
							break;
						}
					}
					if (check) tk = align;
				}

				tk = ap.stime - tk;

				generate_request_function_abstract_tree(sfs,i,tk,ap.d, output);
			}
		}
	}
	timer.end();
	tGenerateRequestFunctionAbstractTree += timer.getTime();
}

void SchedulabilityAnalysis::generate_request_function_abstract_tree(Stateflow** sfs, int i, int stime, int deadline, bool output) {
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		generate_request_function_abstract_tree(sfj,stime,deadline,output);
	}
}

void SchedulabilityAnalysis::generate_request_function_abstract_tree(Stateflow* sf, int stime, int deadline, bool output) {
	//deadline = ceil(1.0*deadline/sf->gcd)*sf->gcd;

	if (sf->mmARFT.find(stime) != sf->mmARFT.end()) {
		if (sf->mmARFT[stime].find(deadline) != sf->mmARFT[stime].end())
			return;
	}

	AbstractRequestFunctionTree arft = AbstractRequestFunctionTree();
	int gcd = Utility::math_gcd(sf->gcd,deadline);
	arft.generate(sf->CRF,stime,deadline,gcd,output);

	sf->mmARFT[stime][deadline] = arft;
}

void SchedulabilityAnalysis::generate_request_function_abstract_tree_for_each_action_instance(Stateflow** sfs, int n, bool output) {
	Timer timer;
	timer.start();
	for (int i=0; i<n-1; i++) {
		Stateflow* sfi = sfs[i];

		generate_request_function_abstract_tree_for_each_action_instance(sfi);
	}
	timer.end();
	tGenerateRequestFunctionAbstractTree += timer.getTime();
}

void SchedulabilityAnalysis::generate_request_function_abstract_tree_for_each_action_instance(Stateflow* sf) {
	for (map<int,map<int,vector<RequestFunction>>>::iterator mmIter = sf->mmCRF.begin(); mmIter != sf->mmCRF.end(); mmIter++) {
		int tk = mmIter->first;
		for (map<int,vector<RequestFunction>>::iterator mIter = mmIter->second.begin(); mIter != mmIter->second.end(); mIter++) {
			int deadline = mIter->first;
			AbstractRequestFunctionTree arft = AbstractRequestFunctionTree();
			int gcd = Utility::math_gcd(sf->gcd,deadline);
			arft.generate(mIter->second,tk,deadline,gcd,false);

			sf->mmARFT[tk][deadline] = arft;
		}
	}
}

void SchedulabilityAnalysis::generate_one_request_function_abstract_tree(Stateflow** sfs, int n, bool output) {
	Timer timer;
	timer.start();
	for (int i=0; i<n-1; i++) {
		Stateflow* sf = sfs[i];
		generate_one_request_function_abstract_tree(sf);
	}
	timer.end();
	tGenerateRequestFunctionAbstractTree += timer.getTime();
}

void SchedulabilityAnalysis::generate_one_request_function_abstract_tree(Stateflow* sf) {
	//deadline = ceil(1.0*deadline/sf->gcd)*sf->gcd;

	sf->arft = AbstractRequestFunctionTree();
	sf->arft.generate2(sf->CRF,0,sf->CRF.at(0).bound,sf->gcd,false);
}

bool SchedulabilityAnalysis::isAbstract(vector<RFNode*> vec_rfnode) {
	for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
		RFNode* rfnode = *iter;
		if (rfnode->isAbstract) return true;
	}
	return false;
}

bool SchedulabilityAnalysis::exact_sched_analysis(Stateflow** sfs, int n, int choice, bool output, ostream& out) {
	Timer timer;
	timer.start();

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		//sf->write_graphviz(cout);
		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
	}

	generate_critical_action_pair(sfs,n);

	if (output) cout << "Generating request functions" << endl;
	//bool success = SchedulabilityAnalysis::generate_request_function(sfs,n,output);
	SchedulabilityAnalysis::generate_request_function_for_each_action_instance(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_critical_request_function(sfs,n,out);

	if (output) cout << "Generating request function abstract tree" <<endl; 
	//SchedulabilityAnalysis::generate_request_function_abstract_tree(sfs,n,output);
	//SchedulabilityAnalysis::generate_one_request_function_abstract_tree(sfs,n,output);
	SchedulabilityAnalysis::generate_request_function_abstract_tree_for_each_action_instance(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_request_function_abstract_tree(sfs,n,out);

	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		for (map<int, vector<ActionPair>>::iterator mIter = sf->mCAP.begin(); mIter != sf->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;

			bool schedulable = exact_sched_analysis(sfs,i,stime,vec_ap);
			if (!schedulable) {
				timer.end();
				tExactStaticOffset += timer.getTime();
				return false;
			}
		}
	}
	timer.end();
	tExactStaticOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::exact_sched_analysis(Stateflow** sfs, int i, int s, vector<ActionPair> vec_ap) {
	for (vector<ActionPair>::iterator iter = vec_ap.begin(); iter != vec_ap.end(); iter++) {
		bool schedulable = exact_sched_analysis(sfs, i, s, *iter);
		if (!schedulable) {
			cout << "UnSchedulable in Stateflow " << i << "\t ActionPair = <" << (*iter).e << "," << (*iter).d << ">" << endl;	
			return false;
		}
	}
	return true;
}


bool SchedulabilityAnalysis::exact_sched_analysis(Stateflow** sfs, int i, int s, ActionPair ap) {
	if (i==0) {
		if (ap.e <= ap.d) return true;
		else return false;
	} else {
		int gcd = 0;
		for (int j=0; j<=i; j++) gcd = Utility::math_gcd(gcd, sfs[j]->gcd);

		vector<AbstractRequestFunctionTree> vec_arft;

		/*
		// one tree
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			vec_arft.push_back(sfj->arft);
		}
		*/
		
		int tk = 0;
		for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
			bool check = true;
			for (int j=0; j<i; j++) {
				Stateflow* sfj = sfs[j];
				if (align % sfj->hyperperiod != 0) {
					check = false;
					break;
				}
			}
			if (check) tk = align;
		}
		tk = ap.stime - tk;
		
		// mutiple trees 
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			map<int, map<int,AbstractRequestFunctionTree>> mmARFT = sfj->mmARFT;
			map<int, map<int,AbstractRequestFunctionTree>>::iterator iter0 = mmARFT.lower_bound(tk);

			if (iter0 == mmARFT.end()) {
				--iter0;
			} else if (tk!=0) 
				iter0 = --iter0;

			if (iter0 == mmARFT.end()) {
				cerr << "Arrive the end of map!=>1" << endl;
				exit(EXIT_FAILURE);
			}
			map<int,AbstractRequestFunctionTree>::iterator iter1 = iter0->second.lower_bound(ap.stime+ap.d);
			if (iter1 == iter0->second.end()) {
				--iter1;
			}
			

			/*
			if (mmARFT.find(tk) == mmARFT.end()) {
				cerr << "Not generate " << tk << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}

			if (mmARFT[tk].find(ap.d) == mmARFT[tk].end()) {
				cerr << "Not generate " << tk << "," << ap.d << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}
			*/
			AbstractRequestFunctionTree arft = iter1->second;
			vec_arft.push_back(arft);
		}
		
		return exact_sched_analysis(vec_arft,ap.e,tk,ap.d,gcd);
	}
}


bool SchedulabilityAnalysis::exact_sched_analysis(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int stime, int deadline, int gcd) {
	vector<vector<RFNode*>> Omega;
	vector<RFNode*> vec_rf;
	for (vector<AbstractRequestFunctionTree>::iterator iter = vec_arft.begin(); iter != vec_arft.end(); iter++) {
		vec_rf.push_back((*iter).root);
	}
	Omega.push_back(vec_rf);

	return exact_sched_analysis(Omega, wcet, stime, deadline, gcd);
}

bool SchedulabilityAnalysis::exact_sched_analysis(vector<vector<RFNode*>> Omega, int wcet, int stime, int deadline, int gcd) {
	for (vector<vector<RFNode*>>::iterator iter = Omega.begin(); iter != Omega.end(); iter++) {
		vector<RFNode*> vec_rfnode = *iter;
		bool schedulable = exact_sched_analysis(vec_rfnode, wcet, stime, deadline, gcd);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::exact_sched_analysis(vector<RFNode*> vec_rfnode, int wcet, int stime, int deadline, int gcd) {
	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);

	// check the schedulability condition
	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		if (rt == INT_MAX) break;
		int start = *iter;
		
		for (int tprim = start+gcd; tprim <= stime+deadline; tprim+=gcd) {
			double sum = 0;
			if (tprim >= stime) sum += wcet;
			for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
				RFNode* rfnode = *iter;
				sum += rfnode->get_mrf(start,tprim);
			}
			if (sum <= tprim-start) {
				rt = max(rt,sum+start-stime);
				break;
			}
			if (tprim == stime+deadline) rt=INT_MAX;
		}
		/*
		if (false && vec_rfnode.size() == 7 && !isAbstract(vec_rfnode)) {
			double sum = wcet;
			for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
				RFNode* rfnode = *iter;
				sum += rfnode->get_mrf(start,77910);
			}
			cout << "77910 => " << sum << endl;
			cout << "Output all the request functions" << endl;

			for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
				RFNode* rfnode = *iter;
				rfnode->rf.output(cout);
			}
		}

		/*
		bool found = false;
		int tprim = start+gcd;
		while (tprim <= stime+deadline) {
			double sum = 0;
			if (tprim >= stime) sum += wcet;
			for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
				RFNode* rfnode = *iter;
				sum += rfnode->get_mrf(start,tprim);
			}
			if (sum <= tprim-start) {
				rt = max(rt,sum+start-stime);
				found = true;
				break;
			}
			int total = ceil(sum)+start;
			if (total == tprim) tprim += 1;
			tprim = total;
		}
		if (!found) rt=INT_MAX;
		*/
	}

	if (rt <= deadline) return true;

	// we need to consider the unschedulable situation
	if (isAbstract(vec_rfnode)) {
		vector<RFNode*> vec_rfnode1;
		vector<RFNode*> vec_rfnode2;

		bool found = false;
		for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
			RFNode* rfnode = *iter;
			if (rfnode->isAbstract && !found) {
				found = true;
				vec_rfnode1.push_back(rfnode->lnode);
				vec_rfnode2.push_back(rfnode->rnode);
				continue;
			}
			vec_rfnode1.push_back(rfnode);
			vec_rfnode2.push_back(rfnode);
		}

		bool schedulable1 = exact_sched_analysis(vec_rfnode1,wcet, stime, deadline, gcd);
		bool schedulable2 = exact_sched_analysis(vec_rfnode2,wcet, stime, deadline, gcd);
		
		return schedulable1 && schedulable2;

	} else 
		return false;
}

int SchedulabilityAnalysis::calculate_remaining_workload(vector<RFNode*> vec_rfnode, int stime, int gcd) {
	int remain = 0;
	int interval = 0;
	int prev = 0;
	for (int t=gcd; t<=stime; t+=gcd) {
		interval += gcd;
		int curr = 0;
		for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
			RFNode* rfnode = *iter;
			if (rfnode->stime != stime) {
				cerr << "calculating remaining workload error" << rfnode->stime << "\t" << stime <<endl;
				exit(EXIT_FAILURE);
			}
			curr += rfnode->get_mrf(0,t);
		}
		remain += curr-prev;
		prev = curr;
		if (remain <= interval) {
			remain = 0;
			interval = 0;
		}
	}

	return remain;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	// release critical request functions
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		// release abstract request function trees
		for (map<int, map<int,AbstractRequestFunctionTree>>::iterator mmIter = sf->mmARFT.begin(); mmIter != sf->mmARFT.end(); mmIter++) {
			map<int,AbstractRequestFunctionTree> temp = mmIter->second;
			for (map<int,AbstractRequestFunctionTree>::iterator mIter = temp.begin(); mIter != temp.end(); mIter++) {
				AbstractRequestFunctionTree& arft = mIter->second;
				arft.root = NULL;
				while(!arft.cur_queue.empty()) arft.cur_queue.pop();
				for (vector<RFNode*>::iterator iter = arft.record.begin(); iter != arft.record.end(); iter++) {
					delete *iter;
				}
				arft.record.clear();
			}
		}

		sf->CRF.clear();
		sf->mmARFT.clear();
	}
	
	Timer timer;
	timer.start();

	bool output = false;

	if (output) cout << "Generating request functions" << endl;
	bool success = SchedulabilityAnalysis::generate_request_function(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_critical_request_function(sfs,n,cout);

	if (output) cout << "Generating request function abstract tree" <<endl; 
	SchedulabilityAnalysis::generate_request_function_abstract_tree(sfs,n,output);
	if (output) SchedulabilityAnalysis::output_request_function_abstract_tree(sfs,n,cout);

	// request functions and abstract rf trees have been generated in function exact_sched_analysis
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];
		for (map<int, vector<ActionPair>>::iterator mIter = sf->mCAP.begin(); mIter != sf->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;

			bool schedulable = exact_analysis_arbitrary_offset(sfs,i,vec_ap);
			if (!schedulable) {
				timer.end();
				//tExactArbitraryOffset += timer.getTime();
				return false;
			}
		}
	}
	timer.end();
	//tExactArbitraryOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(Stateflow** sfs, int i, vector<ActionPair> vec_ap) {
	for (vector<ActionPair>::iterator iter = vec_ap.begin(); iter != vec_ap.end(); iter++) {
		bool schedulable = exact_analysis_arbitrary_offset(sfs, i, *iter);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(Stateflow** sfs, int i, ActionPair ap) {
	if (i==0) {
		if (ap.e <= ap.d) return true;
		else return false;
	} else {
		int gcd = 0;
		for (int j=0; j<=i; j++) gcd = Utility::math_gcd(gcd, sfs[j]->gcd);

		vector<AbstractRequestFunctionTree> vec_arft;
		int tk =ap.stime;

		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			map<int, map<int,AbstractRequestFunctionTree>> mmARFT = sfj->mmARFT;
			if (mmARFT.find(tk) == mmARFT.end()) {
				cerr << "Not generate " << tk << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}

			if (mmARFT[tk].find(ap.d) == mmARFT[tk].end()) {
				cerr << "Not generate " << tk << "," << ap.d << "abstract request function tree for " << ap.toString() <<endl;
				exit(EXIT_FAILURE);
			}
			AbstractRequestFunctionTree arft = mmARFT[tk][ap.d];
			vec_arft.push_back(arft);
		}

		return exact_analysis_arbitrary_offset(vec_arft, ap.e,ap.d,gcd);
	}
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(vector<AbstractRequestFunctionTree> vec_arft, int wcet, int deadline, int gcd) {
	vector<vector<RFNode*>> Omega;
	vector<RFNode*> vec_rf;
	for (vector<AbstractRequestFunctionTree>::iterator iter = vec_arft.begin(); iter != vec_arft.end(); iter++)
		vec_rf.push_back((*iter).root);
	Omega.push_back(vec_rf);
	
	return exact_analysis_arbitrary_offset(Omega, wcet, deadline, gcd);
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(vector<vector<RFNode*>> Omega, int wcet, int deadline, int gcd) {
	for (vector<vector<RFNode*>>::iterator iter = Omega.begin(); iter != Omega.end(); iter++) {
		vector<RFNode*> vec_rfnode = *iter;
		bool schedulable = exact_analysis_arbitrary_offset(vec_rfnode, wcet, deadline, gcd);
		if (!schedulable) return false;
	}
	return true;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset(vector<RFNode*> vec_rfnode, int wcet, int deadline, int gcd) {
	// check the schedulability condition
	double rt = 0;
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (vector<RFNode*>::iterator iter = vec_rfnode.begin(); iter != vec_rfnode.end(); iter++) {
			RFNode* rfnode = *iter;
			sum += rfnode->get_mrf(t);
		}
		if (sum <= t) return true;
	}
	return false;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_simple_digraph();
		digraphs[i] = sf->simple_digraph;
		//digraphs[i]->write_graphviz(cout);
	}
	
	bool schedulable = digraph_exact_analysis(digraphs,n,choice);

	delete[] digraphs;

	timer.end();
	tExactArbitraryOffsetBySimpleDigraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
		//digraphs[i]->write_graphviz(cout);
	}
	
	bool schedulable = digraph_exact_analysis(digraphs,n,choice);

	delete[] digraphs;

	timer.end();
	tExactArbitraryOffsetByPreciseDigraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::exact_analysis_arbitrary_offset_by_precise_digraphs_timeout(Stateflow** sfs, int n, int choice, int timeout) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
		//digraphs[i]->write_graphviz(cout);
	}
	
	// set timeout = 100 seconds
	//prepare_all_digraphs(digraphs,n,choice);
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();
	}

	generate_critical_vertices(digraphs,n);

	// set tf, denoted as the maximum time length for the special digraph
	int tf_max = INT_MIN;
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i];
		for (vector<Node*>::iterator iter = digraph->cnode_vec.begin(); iter != digraph->cnode_vec.end(); iter++) {
			tf_max = max(tf_max,(*iter)->deadline);
		}
	}

	// prepare request functions
	for (int i=0; i<n; i++) {
		Digraph* digraph = digraphs[i]; 
		if (false) cout << "Generate critical request functions..." << endl;
		digraph->generate_critical_request_function(tf_max,false);
		if (false) cout << "Generate abstract request functions tree..." << endl;
		digraph->generate_abstract_request_function_tree();
	}

	bool schedulable = false;
	future<bool> fut = async( SchedulabilityAnalysis::digraph_exact_analysis_timeout,digraphs,n,choice);
	chrono::seconds span(timeout);
	

	if (fut.wait_for(span)==future_status::timeout) {
		nTimeout++;
	}
	else
		schedulable = fut.get();

	delete[] digraphs;

	timer.end();
	tExactArbitraryOffsetByPreciseDigraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset(Stateflow** sfs, int n, int choice, bool periodicity) {
	Timer timer;
	timer.start();
	
	prepare_all_stateflows2(sfs,n,true,periodicity);

	for (int i=0; i<n; i++) {
		bool schedulable = rbf_analysis_static_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tRBFStaticOffset += timer.getTime();	
			return false;
		}
	}
	
	timer.end();
	tRBFStaticOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;

			int tk = 0;
			for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
				bool check = true;
				for (int j=0; j<i; j++) {
					Stateflow* sfj = sfs[j];
					if (align % sfj->hyperperiod != 0) {
						check = false;
						break;
					}
				}
				if (check) tk = align;
			}
			tk = ap.stime - tk;
			bool schedulable = rbf_analysis_static_offset_index(sfs,i,ap.e,tk,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale by rbf analysis with static offset on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;	
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);
	
	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		int start = *iter;
		/*
		for (int tprim = start+gcd; tprim <= stime+deadline; tprim+=gcd) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			for (int j=0; j<i; j++) {
				Stateflow* sfj = sfs[j];
				double rbfsf = sfj->get_rbf(start,tprim);
				temp += rbfsf;
			}
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				break;
			}
			if (tprim == stime+deadline) return false;
		}
		*/
		bool found = false;
		if (start < max(stime-sfs[i]->gcd,0)) continue;

		int tprim = start+gcd;
		while (tprim <= stime+deadline) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			//double temp = wcet;
			for (int j=0; j<i; j++) {
				Stateflow* sfj = sfs[j];
				double rbfsf = sfj->get_rbf(start,tprim);
				temp += rbfsf;
			}
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				found = true;
				break;
			}
			//cout << tprim << endl;
			int total = ceil(1.0*(temp+start)/gcd)*gcd;
			if (total == tprim) tprim += gcd;
			tprim = total;
		}
		if (!found) return false;
	}
	
	/*
	if (rt == 0) {
		cerr << "Error Here!!!!!!!!" << endl;
		exit(EXIT_FAILURE);
	}
	*/

	if (rt <= deadline) return true;
	return false;
}

double SchedulabilityAnalysis::calculate_remaining_rbf_workload(Stateflow** sfs, int i, int stime, int gcd) {
	double remain = 0;
	int interval = 0;
	double prev = 0;
	for (int t=gcd; t<=stime; t+=gcd) {
		interval += gcd;
		double curr = 0;
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			curr += sfj->get_rbf(0,t);
		}
		remain += curr-prev;
		prev = curr;
		if (remain <= interval) {
			remain = 0;
			interval = 0;
		}
	}

	return remain;
}

bool SchedulabilityAnalysis::ibf_analysis_static_offset(Stateflow** sfs, int n, int choice, bool periodicity) {
	Timer timer;
	timer.start();
	
	prepare_all_stateflows2(sfs,n,false,periodicity);

	for (int i=0; i<n; i++) {
		bool schedulable = ibf_analysis_static_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tIBFStaticOffset += timer.getTime();
			return false;
		}
	}
	
	timer.end();
	tIBFStaticOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_static_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;

			int tk = 0;
			for (int align = 0; align <=ap.stime; align += sfs[i]->hyperperiod) {
				bool check = true;
				for (int j=0; j<i; j++) {
					Stateflow* sfj = sfs[j];
					if (align % sfj->hyperperiod != 0) {
						check = false;
						break;
					}
				}
				if (check) tk = align;
			}
			tk = ap.stime - tk;
			bool schedulable = ibf_analysis_static_offset_index(sfs,i,ap.e,tk,ap.d);
			if (!schedulable) return false;
		}
	}
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_static_offset_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	int t_gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
		t_gcd = Utility::math_gcd(t_gcd,sfj->t_gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	// generate the start times for all the busy periods
	set<int> startSet;
	for (int start = 0; start <=stime; start+=gcd) 
		startSet.insert(start);

	/*
	if (i==7) {
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			cout << "sf_" << j << ".ibf[0,77910)=" << sfj->get_ibf(0,77910) << endl;
		}
	}
	*/
	
	double rt = 0;
	for (set<int>::iterator iter = startSet.begin(); iter != startSet.end(); iter++) {
		int start = *iter;
		bool found = false;
		if (start < max(stime-sfs[i]->gcd,0)) continue;
		int tprim = start+t_gcd;
		while (tprim <= stime+deadline) {
			double temp = 0;
			if (tprim >= stime) temp += wcet;
			//double temp = wcet;
			for (int j=0; j<i; j++) {
				Stateflow* sfj = sfs[j];
				double ibfsf = sfj->get_ibf(start,tprim);
				temp += ibfsf;
			}
			if (temp <= tprim-start) {
				rt = max(rt,temp+start-stime);
				found = true;
				break;
			}
			int total = ceil(temp)+start;
			if (total == tprim) tprim += t_gcd;
			tprim = total;
		}
		if (!found) return false;
	}

	if (rt <= deadline) return true;
	return false;
}

double SchedulabilityAnalysis::calculate_remaining_ibf_workload(Stateflow** sfs, int i, int stime, int gcd) {
	return calculate_remaining_rbf_workload(sfs,i,stime,gcd);
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();

	prepare_all_stateflows2(sfs,n,true,false);

	for (int i=0; i<n; i++) {
		bool schedulable = rbf_analysis_arbitrary_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tRBFArbitraryOffset += timer.getTime();
			return false;
		}
	}
	
	timer.end();
	tRBFArbitraryOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			bool schedulable = rbf_analysis_arbitrary_offset_index(sfs,i,ap.e,ap.d);
			if (!schedulable) {
				cout << "Unschedulable by RBF analysis with arbitrary offset for Stateflow " << i 
					<< "\t Action Pair=>" << ap.toString() <<endl;	
				return false;
			}
		}
	}
	
	return true;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	// check the schedulability condition
	for (int t=0; t<=deadline; t+=gcd) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			sum += sfj->get_rbf(t);
		}
		if (sum <= t) return true;
	}

	return false;
}



bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice, bool linearization) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_simple_digraph();
		digraphs[i] = sf->simple_digraph;
		//digraphs[i]->write_graphviz(cout);
	}
	
	bool schedulable = digraph_rbf_analysis(digraphs,n,choice,linearization);

	delete[] digraphs;

	timer.end();
	if (linearization)
		tRBFArbitraryOffsetBySimpleDigraphWithLinearization += timer.getTime();
	else
		tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_simple_digraphs_DP(Stateflow** sfs, int n, int choice, bool periodicityProperty) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_simple_digraph();
		digraphs[i] = sf->simple_digraph;
		//digraphs[i]->write_graphviz(cout);
	}
	
	bool schedulable = digraph_rbf_analysis_dynamic_programming(digraphs,n,choice,periodicityProperty,false);

	delete[] digraphs;

	timer.end();

	tRBFArbitraryOffsetBySimpleDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice, bool linearization) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
	}

	bool schedulable = digraph_rbf_analysis(digraphs,n,choice,linearization);

	delete[] digraphs;

	timer.end();
	if (linearization)
		tRBFArbitraryOffsetByPreciseDigraphWithLinearization += timer.getTime();
	else
		tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_precise_digraphs_DP(Stateflow** sfs, int n, int choice,bool periodicityProperty)  {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
	}

	bool schedulable = digraph_rbf_analysis_dynamic_programming(digraphs,n,choice,periodicityProperty,false);

	delete[] digraphs;

	timer.end();
	
	tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_strongly_connected_precise_digraphs_DP(Stateflow** sfs, int n, int choice,bool periodicityProperty)  {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_strongly_connected_precise_digraph();
		digraphs[i] = sf->sc_precise_digraph;

		sf->generate_precise_digraph();
		
		int num = 0;
		for (vector<Node*>::iterator iter = sf->precise_digraph->node_vec.begin(); iter != sf->precise_digraph->node_vec.end(); iter++) {
			if ((*iter)->in.empty())
				num++;
		}
		cout << "Precise: nNode=" << sf->precise_digraph->node_vec.size() << ", nEdge=" << sf->precise_digraph->edge_vec.size() << ", nIsolate=" << num << endl;
		cout << "SC_Precise: nNode=" << sf->sc_precise_digraph->node_vec.size() << ", nEdge=" << sf->sc_precise_digraph->edge_vec.size() << endl;
		//sf->precise_digraph->write_graphviz(cout);
		//digraphs[i]->write_graphviz(cout);
	}

	bool schedulable = digraph_rbf_analysis_dynamic_programming(digraphs,n,choice,periodicityProperty,false);

	delete[] digraphs;

	timer.end();
	
	tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::rbf_analysis_arbitrary_offset_by_reachability_digraphs_DP(Stateflow** sfs, int n, int choice,bool periodicityProperty)  {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_reachability_digraph();
		digraphs[i] = sf->reachability_digraph;
	}

	bool schedulable = reachability_digraph_rbf_analysis_dynamic_programming(digraphs,n,choice);

	delete[] digraphs;

	timer.end();
	
	tRBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();

	prepare_all_stateflows2(sfs,n,false,false);

	for (int i=0; i<n; i++) {
		Stateflow* sfi = sfs[i];
		sfi->calculate_deadlines();
		bool schedulable = ibf_analysis_arbitrary_offset_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tIBFArbitraryOffset += timer.getTime();
			return false;
		}
	}
	
	timer.end();
	tIBFArbitraryOffset += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			bool schedulable = ibf_analysis_arbitrary_offset_index(sfs,i,ap.e,ap.d);
			if (!schedulable) {
				cout << "Unschedulable by IBF analysis with arbitrary offset for Stateflow " << i 
					<< "\t Action Pair=>" << ap.toString() <<endl;	
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_index(Stateflow** sfs, int i, int wcet, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->t_gcd);
	}
	
	// check the schedulability condition
	int t = 0;
	while (t<=deadline) {
		double sum = wcet;
		for (int j=0; j<i; j++) {
			Stateflow* sfj = sfs[j];
			sum += sfj->get_ibf(t);
		}
		if (sum <= t) return true;
		int total = ceil(sum);
		if (total==t) t+=gcd;
		t = total;
	}

	return false;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs(Stateflow** sfs, int n, int choice, bool linearization) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_simple_digraph();
		digraphs[i] = sf->simple_digraph;
	}

	bool schedulable = digraph_ibf_analysis(digraphs,n,choice,linearization);

	delete[] digraphs;

	timer.end();
	if (linearization)
		tIBFArbitraryOffsetBySimpleDigraphWithLinearization += timer.getTime();
	else
		tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_simple_digraphs_DP(Stateflow** sfs, int n, int choice,bool periodicityProperty) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_simple_digraph();
		digraphs[i] = sf->simple_digraph;
	}

	bool schedulable = digraph_ibf_analysis_dynamic_programming(digraphs,n,choice,periodicityProperty,false);

	delete[] digraphs;

	timer.end();
	
	tIBFArbitraryOffsetBySimpleDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs(Stateflow** sfs, int n, int choice, bool linearization) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
	}

	bool schedulable = digraph_ibf_analysis(digraphs,n,choice,linearization);

	delete[] digraphs;

	timer.end();
	if (linearization)
		tIBFArbitraryOffsetByPreciseDigraphWithLinearization += timer.getTime();
	else
		tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_precise_digraphs_DP(Stateflow** sfs, int n, int choice,bool periodicityProperty) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_precise_digraph();
		digraphs[i] = sf->precise_digraph;
	}

	bool schedulable = digraph_ibf_analysis_dynamic_programming(digraphs,n,choice,periodicityProperty,false);

	delete[] digraphs;

	timer.end();
	
	tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_strongly_connected_precise_digraphs_DP(Stateflow** sfs, int n, int choice,bool periodicityProperty) {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_strongly_connected_precise_digraph();
		digraphs[i] = sf->sc_precise_digraph;
	}

	bool schedulable = digraph_ibf_analysis_dynamic_programming(digraphs,n,choice,periodicityProperty,false);

	delete[] digraphs;

	timer.end();
	
	tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::ibf_analysis_arbitrary_offset_by_reachability_digraphs_DP(Stateflow** sfs, int n, int choice,bool periodicityProperty)  {
	Timer timer;
	timer.start();
	Digraph** digraphs = new Digraph*[n];
	// Generate all the simple digraphs for stateflows
	for (int i=0; i<n; i++) {
		Stateflow* sf = sfs[i];

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();

		sf->generate_reachability_digraph();
		digraphs[i] = sf->reachability_digraph;
	}

	bool schedulable = reachability_digraph_ibf_analysis_dynamic_programming(digraphs,n,choice);

	delete[] digraphs;

	timer.end();
	
	tIBFArbitraryOffsetByPreciseDigraphWithSearchGraph += timer.getTime();
	return schedulable;
}

bool SchedulabilityAnalysis::lu_rbf_sched_analysis(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();

	prepare_all_stateflows(sfs,n,choice);

	for (int i=0; i<n; i++) {
		bool schedulable = lu_rbf_sched_analysis_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tLinearUpperRBF += timer.getTime();
			return false;
		}
	}
	timer.end();
	tLinearUpperRBF += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::lu_rbf_sched_analysis_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			
			bool schedulable = lu_rbf_sched_analysis_index(sfs,i,ap.e,ap.stime,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::lu_rbf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	double tUtil = 0;
	double tCrbf = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCrbf += sfj->crbf;
		if (tUtil >= 1) return false;
	}

	double rt = (tCrbf+wcet)/(1.0-tUtil);
	if (rt <= deadline) return true;
	return false;
}

bool SchedulabilityAnalysis::lu_ibf_sched_analysis(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();

	prepare_all_stateflows(sfs,n,choice);

	for (int i=0; i<n; i++) {
		bool schedulable = lu_ibf_sched_analysis_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tLinearUpperIBF += timer.getTime();
			return false;
		}
	}
	timer.end();
	tLinearUpperIBF += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::lu_ibf_sched_analysis_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			
			bool schedulable = lu_ibf_sched_analysis_index(sfs,i,ap.e,ap.stime,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::lu_ibf_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	double tUtil = 0;
	double tCibf = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCibf += sfj->cibf;
		if (tUtil >= 1) return false;
	}

	double rt = (tCibf+wcet)/(1.0-tUtil);
	if (rt <= deadline) return true;
	return false;
}

bool SchedulabilityAnalysis::lu_csum_sched_analysis(Stateflow** sfs, int n, int choice) {
	Timer timer;
	timer.start();

	prepare_all_stateflows(sfs,n,choice);

	for (int i=0; i<n; i++) {
		bool schedulable = lu_csum_sched_analysis_index(sfs,i);
		if (!schedulable) {
			timer.end();
			tLinearUpperCSUM += timer.getTime();
			return false;
		}
	}
	timer.end();
	tLinearUpperCSUM += timer.getTime();
	return true;
}

bool SchedulabilityAnalysis::lu_csum_sched_analysis_index(Stateflow** sfs, int i) {
	Stateflow* sfi = sfs[i];
	for (map<int,vector<ActionPair>>::iterator mIter = sfi->mCAP.begin(); mIter != sfi->mCAP.end(); mIter++) {
		int stime = mIter->first;
		vector<ActionPair> vec_ap = mIter->second;
		for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
			ActionPair ap = *apIter;
			
			bool schedulable = lu_csum_sched_analysis_index(sfs,i,ap.e,ap.stime,ap.d);
			if (!schedulable) {
				cout << "UnSchedulbale on Stateflow " << i << "\t Action Pair=>" << ap.toString() <<endl;
				return false;
			}
		}
	}
	return true;
}

bool SchedulabilityAnalysis::lu_csum_sched_analysis_index(Stateflow** sfs, int i, int wcet, int stime, int deadline) {
	if (i==0) {
		if (wcet <= deadline) return true;
		else return false;
	}

	int gcd = 0;
	for (int j=0; j<=i; j++) {
		Stateflow* sfj = sfs[j];
		gcd = Utility::math_gcd(gcd,sfj->gcd);
	}

	if (stime%gcd != 0) {
		cerr << "stime=" << stime << " cannot be divided by gcd=" << gcd <<endl;
		exit(EXIT_FAILURE);
	}

	double tUtil = 0;
	double tCSum = 0;
	for (int j=0; j<i; j++) {
		Stateflow* sfj = sfs[j];
		tUtil += sfj->lfac;
		tCSum += 2.0*sfj->csum;
		if (tUtil >= 1) return false;
	}

	double rt = (tCSum+wcet)/(1.0-tUtil);
	if (rt <= deadline) return true;
	return false;
}

void SchedulabilityAnalysis::output_critical_action_pair(Stateflow** sfs, int n, ostream& out) {
	for (int i=0; i<n; i++) {
		out << "Stateflow "<<i<< " "": Show critical action pairs:"<<endl;
		Stateflow* sf = sfs[i];

		for (map<int,vector<ActionPair>>::iterator mIter = sf->mCAP.begin(); mIter != sf->mCAP.end(); mIter++) {
			int stime = mIter->first;
			vector<ActionPair> vec_ap = mIter->second;

			out << "\t" << "Trigger time = "<<stime<<endl;
			for (vector<ActionPair>::iterator apIter = vec_ap.begin(); apIter != vec_ap.end(); apIter++) {
				ActionPair ap = *apIter;
				out << "\t\t" << "<" << ap.priority << "," << ap.stime << "," << ap.e << "," << ap.d << ">" <<endl; 
			}
		}
		out <<endl;
	}
}

void SchedulabilityAnalysis::output_critical_request_function(Stateflow** sfs, int n, ostream& out) {
	for (int i=0; i<n; i++) {
		out << "Stateflow " <<i<<": Show critical request functions:"<<endl;
		Stateflow* sf = sfs[i];

		output_critical_request_function(sf,out);
	}
}

void SchedulabilityAnalysis::output_critical_request_function(Stateflow* sf, ostream& out) {
	out << "Number=" << sf->CRF.size() <<endl;
	output_critical_request_function(sf->CRF,out);
	out <<endl;
}

void SchedulabilityAnalysis::output_critical_request_function(vector<RequestFunction> vec_rf, ostream& out) {
	for (vector<RequestFunction>::iterator rfIter = vec_rf.begin(); rfIter != vec_rf.end(); rfIter++) {
		RequestFunction rf = *rfIter;
		rf.output(out);
	}

}

void SchedulabilityAnalysis::output_request_function_abstract_tree(Stateflow** sfs, int n, ostream& out) {
	for (int i=0; i<n-1; i++) {
		out << "Stateflow " <<i<<": Show request functions abstract tree:"<<endl;
		Stateflow* sf = sfs[i];

		output_request_function_abstract_tree(sf,out);
	}
}

void SchedulabilityAnalysis::output_request_function_abstract_tree(Stateflow* sf, ostream& out) {
	for (map<int, map<int,AbstractRequestFunctionTree>>::iterator mIter = sf->mmARFT.begin(); mIter != sf->mmARFT.end(); mIter++) {
		int stime = mIter->first;
		map<int,AbstractRequestFunctionTree> mARFT = mIter->second;
		for (map<int,AbstractRequestFunctionTree>::iterator mIter2 = mARFT.begin(); mIter2 != mARFT.end(); mIter2++) {
			int deadline = mIter2->first;
			AbstractRequestFunctionTree arft = mIter2->second;
			out << "Within [" << stime << "," << deadline <<")" <<endl;
			output_request_function_abstract_tree(arft,out);
			out << endl;
		}
	}
}

void SchedulabilityAnalysis::output_request_function_abstract_tree(AbstractRequestFunctionTree arft, ostream& out) {
	arft.root->output(out);
}