#ifndef DIGRAPHREQUESTFUNCTION_H_
#define DIGRAPHREQUESTFUNCTION_H_

#include "NodeAndEdge.h"

class DigraphRequestFunction {
public:
	int priority;
	int ub; // upper bound on response time
	int gcd;

	int lastValue; // the last value
	int lastTime; // the trigger time of the end node
	Node* lastNode; // the last edge
	bool isCounted; // whether the last transition has been counted 

	list<Node*> path; // list of transitions
	map<int,int> rf; // rf(t)

	DigraphRequestFunction() {}
	~DigraphRequestFunction() {
		//path.clear();
		//rf.clear();
	}

	DigraphRequestFunction(int _priority, int _ub, int _gcd) {
		priority = _priority;
		ub = _ub;
		gcd = _gcd;

		lastValue = 0;
		lastNode = NULL;
	}

	DigraphRequestFunction(const DigraphRequestFunction& _rf) {
		priority = _rf.priority;
		ub = _rf.ub;
		gcd = _rf.gcd;

		lastValue = _rf.lastValue;
		lastNode = _rf.lastNode;
		isCounted = _rf.isCounted;

		path = _rf.path;
		rf = _rf.rf;
	}

	void add_node (int curTime, Node* node) {
		if (node == NULL) return;

		path.push_back(node);

		rf[curTime] = lastValue;
		lastValue += node->wcet;
		lastTime = curTime;
		lastNode = node;
	}

	int request_function(int t) {
		if (t > lastTime) return lastValue;

		map<int,int>::iterator iter = rf.lower_bound(t);

		return iter->second;
	}

	const bool operator>(DigraphRequestFunction& rf) {
		// the compared rf should belong to the same digraph
		if (priority != rf.priority) return false;
		if (gcd != rf.gcd) return false;

		if (lastValue < rf.lastValue) return false;

		int finish = max(lastTime, rf.lastTime);

		for (int t = 0; t<=finish; t+=gcd)
			if (request_function(t) < rf.request_function(t)) return false;

		return true;
	}

	void output(ostream& out) {
		out << "========================================" << endl;
		out << "RF-PATH" << endl;
		for (list<Node*>::iterator iter = path.begin(); iter != path.end(); iter++) {
			out << (*iter)->name;
			if (iter != (--path.end())) out << "->";
		}
		out << endl;
		
		out << "RF-Value" << endl;
		for (map<int,int>::iterator iter = rf.begin(); iter != rf.end(); iter++) {
			out << "rf(" << iter->first << ")=" << iter->second << endl;
		}
		out << "========================================" << endl;
	}
};

#endif