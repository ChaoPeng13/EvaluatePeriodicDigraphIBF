#include "RandomGenerator.h"
#include "GraphAlgorithms.h"
#include "FileReader.h"

// Initialize Digraph Generating Parameters
int	RandomGenerator::DIGRAPH_SCALE = 1000;
int RandomGenerator::DIGRAPH_PERIOD[3][2] = {{2,4},{6,12},{5,10}};
int RandomGenerator::DIGRAPH_SCALE_FACTOR[5] = {1,2,4,5,10};
int RandomGenerator::DIGRAPH_SCALE_FACTOR2[4] = {1,10,100,1000};

// Initialize Stateflow Generating Parameters
int RandomGenerator::STATEFLOW_SCALE = 1000;
int RandomGenerator::STATEFLOW_PERIOD[6] = {5, 10, 20, 25, 50, 100};
int RandomGenerator::STATEFLOW_PERIOD2[5] = {10, 20, 25, 50, 100};
int RandomGenerator::STATEFLOW_PERIOD3[4] = {10, 20, 25, 50};
int RandomGenerator::STATEFLOW_PERIOD4[5] = {5, 10, 20, 25, 50};
int RandomGenerator::STATEFLOW_PERIOD5[4] = {1, 2, 5, 10};
int RandomGenerator::STATEFLOW_PERIOD6[8] = {1,2,5,10,20,25,50,100};
int RandomGenerator::STATEFLOW_PERIOD7[3] = {10, 20, 40};
int RandomGenerator::STATEFLOW_PERIOD8[10] = {1,2,5,10,20,25,50,100,200,500};
int RandomGenerator::STATEFLOW_PERIOD9[10] = {1,2,5,10,20,25,50,100,200,250};
int RandomGenerator::STATEFLOW_PERIOD10[6] = {10, 20, 25, 40, 50, 100};
int RandomGenerator::STATEFLOW_PERIOD11[4] = {10, 20, 40, 80};
int RandomGenerator::STATEFLOW_PERIOD12[10] = {1,2,5,10,20,25,50,100,200,1000};
int RandomGenerator::STATEFLOW_PERIOD13[9] = {1,2,5,10,20,25,50,100,1000};
int RandomGenerator::STATEFLOW_PERIOD14[9] = {1,2,5,10,20,50,100,1000,2000};
int RandomGenerator::STATEFLOW_PERIOD15[5] = {10, 20, 25,50, 100};
int RandomGenerator::STATEFLOW_PERIOD16[4] = {10, 20, 50, 100};
int RandomGenerator::STATEFLOW_PERIOD17[10] = {1,2,5,10,20,25,50,100,1000,2000};
//int RandomGenerator::STATEFLOW_FACTOR[6] = {1, 2, 4, 5, 10, 20};
int RandomGenerator::STATEFLOW_BASE[8] = {1100, 1300,2300, 3100, 3700, 3000, 5000, 7000};
//int RandomGenerator::STATEFLOW_BASE[5] = {10000, 11000, 13000, 15000, 17000};
//int RandomGenerator::STATEFLOW_BASE[8] = {110, 130,230, 310, 370, 300, 500, 700};
int RandomGenerator::STATEFLOW_FACTOR[5] = {1, 2, 4, 5, 10};
int RandomGenerator::STATEFLOW_FACTOR2[4] = {1, 2, 4, 5};
int RandomGenerator::STATEFLOW_FACTOR3[3] = {1, 2, 4};
int RandomGenerator::STATEFLOW_FACTOR4[2] = {1, 2};

Digraph* RandomGenerator::generate_one_digraph(int index, int scale, int nNode, int nEdge) {
	Digraph* digraph = new Digraph(index, scale);

	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	connected.push_back(nodes[rand()%nNode]);
	for (int i=0; i<nNode; i++) {
		Node* newNode = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == newNode) {
				isContained = true;
				break;
			}
		}

		if (isContained) continue;

		Node* connectedNode = connected.at(rand()%connected.size());
		Edge* edge = NULL;

		int new2Connected = rand()%2;

		if (new2Connected) 
			edge = new Edge(newNode, connectedNode);
		else 
			edge = new Edge(connectedNode, newNode);

		digraph->add_edge(edge);
		connected.push_back(newNode);
	}

	for (int i= nNode-1; i<nEdge; i++) {
		Node* src = nodes[rand()%nNode];
		Node* snk = nodes[rand()%nNode];

		Edge* edge = new Edge(src,snk);
		digraph->add_edge(edge);
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		for (int i=0; i<nEdge; i++) {
			int reverse = rand()%2;

			if(reverse) {
				Edge* edge = digraph->edge_vec.at(i);

				// reverse the direction of the edge
				edge->src_node->in.push_back(edge);
				edge->src_node->out.remove(edge);

				edge->snk_node->out.push_back(edge);
				edge->snk_node->in.remove(edge);

				Node* temp = edge->src_node;
				edge->src_node = edge->snk_node;
				edge->snk_node = temp;
			}
		}
	}

	return digraph;
}

Digraph* RandomGenerator::generate_one_digraph(int index, int scale, int nNode) {
	Digraph* digraph = new Digraph(index, scale);

	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	delete[] nodes;

	return digraph;
}

Digraph RandomGenerator::generate_one_digraph_without_pointer(int index, int scale, int nNode) {
	Digraph digraph(index, scale);

	vector<Node> node_vec;
	vector<Edge> edge_vec;

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node node(name, i, scale);
		node_vec.push_back(node);
	}

	set<int> connected;
	for (int i=0; i<nNode; i++) {
		Node &node = node_vec[i];

		if (connected.find(i) == connected.end()) connected.insert(i);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node &newNode = node_vec[start+j];

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) {
				Edge edge(&node,&newNode);
				
				node.out.push_back(&edge);
				newNode.in.push_back(&edge);

				edge_vec.push_back(edge);
			}
			else {
				Edge edge(&newNode,&node);

				node.in.push_back(&edge);
				newNode.out.push_back(&edge);

				edge_vec.push_back(edge);
			}
		}
	}

	Node& src = node_vec[0];
	for (int i=0; i<nNode; i++) {
		Node &node = node_vec[i];

		// adding backedge from leaf node to src
		if (node.out.size() == 0) {
			Edge edge(&node,&src);
			edge_vec.push_back(&edge);
		}
	}

	while(true) {
		digraph.update(node_vec,edge_vec);
		Digraph* pDigraph = &digraph;
		if (GraphAlgorithms::isCyclic(pDigraph)) break;

		Edge& edge = edge_vec.at(rand()%edge_vec.size());

		// reverse the direction of the edge
		edge.src_node->in.push_back(&edge);
		edge.src_node->out.remove(&edge);

		edge.snk_node->out.push_back(&edge);
		edge.snk_node->in.remove(&edge);

		Node* temp = edge.src_node;
		edge.src_node = edge.snk_node;
		edge.snk_node = temp;
	}

	return digraph;
}

Digraph* RandomGenerator::generate_one_digraph(int index, int scale, int nNode, double util) {
	 Digraph* digraph = RandomGenerator::generate_one_digraph(0,RandomGenerator::DIGRAPH_SCALE,nNode);
			
	int base =  calculate_base();

	// generate random edges with random separation time
	for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
		Edge* edge = *iter;

		int separationTime = calculate_separation_time(base) * 100;
		edge->set_separation_time(separationTime);
	}

	// generate random nodes with random deadline and wcet
	// set constrained deadlines
	for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
		Node* node = *iter;

		int minSeparationTime = INT_MAX;
		for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) 
			minSeparationTime = min(minSeparationTime, (*e)->separationTime);

		int deadline = minSeparationTime;
		node->set_deadline(deadline);

		//int wcet = rand()%(deadline/2)+1;
		//int wcet = 1.0*(rand()%1000)*deadline;
		int wcet = util*deadline;
		if (wcet == 0) wcet = 1;
		node->set_wcet(wcet);
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	digraph->generate_strongly_connected_components();
	digraph->check_strongly_connected();
	digraph->calculate_period_gcd();
	digraph->calculate_linear_factor();

	cout<<"linear factor: expected="<<util<<"\tfact="<<digraph->linear_factor<<endl;
	// step 2: rescale wcets
	int times = 0;
	while (abs(util-digraph->linear_factor) > 0.001) {
		double factor = util/digraph->linear_factor;
		
		// release sccs
		for (vector<Digraph*>::iterator iter = digraph->sccs.begin(); iter != digraph->sccs.end(); iter++) {
			(*iter)->node_vec.clear();
			(*iter)->edge_vec.clear();
			delete *iter;
			*iter = NULL;
		}
		digraph->sccs.clear();

		digraph->scale_wcet(factor);
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		digraph->calculate_linear_factor();

		cout << "Run-"<< ++times <<"\tlinear factor: expected="<<util<<"\tfact="<<digraph->linear_factor<<endl;
	}

	return digraph;
}

Digraph* RandomGenerator::generate_one_digraph2(int index, int scale, int nNode) {
	Digraph* digraph = new Digraph(index, scale);

	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = {10,20}
	// wcets = {1,2,5}
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		//edge->separationTime = 10*(1+rand()%2);
		edge->separationTime = 10+rand()%11;
	}
	
	int wcets[3] = {1,2,5};
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		node->wcet = wcets[rand()%3];
	}

	return digraph;
}

Digraph* RandomGenerator::generate_one_digraph3(int index, int scale) {
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [5,10]
	int nNode = 5+rand()%6;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = [100,300]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		//edge->separationTime = 10*(1+rand()%2);
		//edge->separationTime = 100+rand()%201;
		edge->separationTime = 100+rand()%101;
		//edge->separationTime = 10+rand()%91;
		//edge->separationTime = 200+rand()%801;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// wcets = [0,0.07]*deadline
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		double r = 0.07*(rand()%100+1)/100;
		node->wcet = ceil(r*node->deadline);
	}

	return digraph;
}

Digraph* RandomGenerator::generate_one_digraph4(int index, int scale, double util) {
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [5,10]
	int nNode = 5+rand()%6;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = [100,300]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		//edge->separationTime = 10*(1+rand()%2);
		//edge->separationTime = 100+rand()%201;
		//edge->separationTime = (100+rand()%101)*;
		//edge->separationTime = (10+rand()%91)*scale;
		//edge->separationTime = 200+rand()%801;
		edge->separationTime = (50+rand()%101)*scale;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// wcets = util*deadline
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int wcet = floor(util*node->deadline);
		if (wcet == 0) wcet = 1;
		node->wcet = wcet;
	}

	digraph->generate_strongly_connected_components();
	digraph->calculate_linear_factor();

	cout << index << "\tExpectedUtil=" << util << "\tRealUtil=" << digraph->linear_factor << endl;
	/*
	int run =0;
	while (abs(digraph->linear_factor-util)>0.01) {
		double factor = util/digraph->linear_factor;
		digraph->scale_wcet(factor);

		// release sccs
		for (vector<Digraph*>::iterator iter = digraph->sccs.begin(); iter != digraph->sccs.end(); iter++) {
			(*iter)->node_vec.clear();
			(*iter)->edge_vec.clear();
			delete *iter;
			*iter = NULL;
		}
		digraph->sccs.clear();

		digraph->generate_strongly_connected_components();
		digraph->calculate_linear_factor();
		cout << index << "\tRun:" << run++ << "\tExpectedUtil=" << util << "\tRealUtil=" << digraph->linear_factor << endl;
	}
	*/
	return digraph;
}

Digraph* RandomGenerator::generate_DIGRAPH_I(int index, double util) {
	int scale = 1;
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [5,10]
	int nNode = 5+rand()%6;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = [20,100]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		edge->separationTime = (20+rand()%81)*scale;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// wcets = util*deadline
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int wcet = floor(util*node->deadline);
		if (wcet == 0) wcet = 1;
		node->wcet = wcet;
	}
	
	delete[] nodes; 
	return digraph;
}

Digraph* RandomGenerator::generate_DIGRAPH_II(int index, double util) {
	int SCALE_SET[3] = {10,100,1000}; 
	int scale = SCALE_SET[rand()%3];

	Digraph* digraph = new Digraph(index, scale);

	// vertices = [5,10]
	int nNode = 5+rand()%6;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = [100,200]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		edge->separationTime = (100+rand()%101)*scale;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// wcets = util*deadline
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int wcet = floor(util*node->deadline);
		if (wcet == 0) wcet = 1;
		node->wcet = wcet;
	}
	
	delete[] nodes;
	return digraph;
}

Digraph* RandomGenerator::generate_DIGRAPH_III(int index) {
	int scale = 1;
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [1,10]
	int nNode = 1+rand()%10;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = [50,100]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		edge->separationTime = (50+rand()%51)*scale;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// wcets = [1,8]
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		node->wcet = 1+rand()%8;
	}
	
	delete[] nodes; 
	return digraph;
}

Digraph* RandomGenerator::generate_DIGRAPH_IV(int index,double util) {
	int scale = 1;
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [5,10]
	int nNode = 5+rand()%6;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = [50,100]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		edge->separationTime = (50+rand()%51)*scale;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);
		node->wcet = util*deadline;

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	
	delete[] nodes; 
	return digraph;
}

Digraph* RandomGenerator::generate_DIGRAPH_V(int index,double util) {
	int scale = 100;
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [1,10]
	int nNode = 1+rand()%10;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// periods = [200,500]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		edge->separationTime = (200+rand()%301)*scale;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);
		node->wcet = ceil(util*deadline);

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}

	// scale the wcets
	digraph->generate_strongly_connected_components();
	digraph->check_strongly_connected();
	digraph->calculate_period_gcd();
	digraph->calculate_all_gcd();
	digraph->calculate_linear_factor();

	double factor = util/digraph->linear_factor;
	digraph->scale_wcet(factor);

	// release the sccs
	for (vector<Digraph*>::iterator iter = digraph->sccs.begin(); iter != digraph->sccs.end(); iter++) {
		(*iter)->node_vec.clear();
		(*iter)->edge_vec.clear();
		delete *iter;
		*iter = NULL;
	}
	digraph->sccs.clear();
	
	delete[] nodes; 
	return digraph;
}

Digraph* RandomGenerator::generate_DIGRAPH_VI(int index) {
	int scale = 1;
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [7,15]
	int nNode = 7+rand()%9;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	/// two types of verticies:
	/// Type-1: wcet = [1,6] and period = [20,100]
	/// Type-2: wcet = [7,9] and period = [101,300]
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int r = rand()%2;

		if (r==0) {
			node->wcet = 1+rand()%6;
		}
		else {
			node->wcet = 7+rand()%3;
		}

		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;

			if (r==0) {
				edge->separationTime = (20+rand()%81)*scale;
			}
			else {
				edge->separationTime = (101+rand()%200)*scale;
			}

			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);

		if (node->out.empty()) {
			if (r==0) 
				node->set_deadline(100);
			else node->set_deadline(300);
		}
	}

	delete[] nodes; 
	return digraph;
}

Digraph* RandomGenerator::generate_DIGRAPH_VII(int index, int periodType) {
	int scale = 1;
	Digraph* digraph = new Digraph(index, scale);

	// vertices = [5,7]
	int nNode = 5+rand()%3;
	Node** nodes = new Node*[nNode];

	for (int i=0; i<nNode; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		Node* node = new Node(name, i, scale);
		nodes[i] = node;
		digraph->add_node(node);
	}

	vector<Node*> connected;
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		bool isContained = false;
		for (vector<Node*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == node) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(node);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int numEdge = 0;
		if(prob <= 0.4)      numEdge = 1;
		else if(prob <= 0.8) numEdge = 2;
		else if(prob <= 0.9) numEdge = 3;
		else if(prob <= 1.0) numEdge = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		numEdge = min(numEdge, nNode-start);
		for (int j=0; j<numEdge; j++) {
			Node* newNode = nodes[start+j];
			Edge* edge = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				edge = new Edge(node,newNode);
			else
				edge = new Edge(newNode,node);

			digraph->add_edge(edge);
		}
	}

	Node* src = nodes[0];
	for (int i=0; i<nNode; i++) {
		Node* node = nodes[i];

		// adding backedge from leaf node to src
		if (node->out.size() == 0) {
			Edge* edge = new Edge(node,src);
			digraph->add_edge(edge);
		}
	}

	while(!GraphAlgorithms::isCyclic(digraph)) {
		Edge* edge = digraph->edge_vec.at(rand()%digraph->edge_vec.size());

		// reverse the direction of the edge
		edge->src_node->in.push_back(edge);
		edge->src_node->out.remove(edge);

		edge->snk_node->out.push_back(edge);
		edge->snk_node->in.remove(edge);

		Node* temp = edge->src_node;
		edge->src_node = edge->snk_node;
		edge->snk_node = temp;
	}

	// set up the wcets and periods
	// Type-1: period = [20,40]
	// Type-2: and period = [40,60]
	// Type-3: and period = [60,80]
	for (vector<Edge*>::iterator eIter = digraph->edge_vec.begin(); eIter != digraph->edge_vec.end(); eIter++) {
		Edge* edge = *eIter;
		int separationTime;
		if (periodType == 1) 
			separationTime = (20+rand()%21)*scale;
		else if (periodType == 2)
			separationTime = (40+rand()%21)*scale;
		else if (periodType == 3)
			separationTime = (60+rand()%21)*scale;
		else {
			cerr << "Error the period type!\t" << periodType << endl;
			exit(EXIT_FAILURE);
		}

		edge->separationTime = separationTime;
	}
	
	// deadline = the minimum inter-release time of all the outgoing edges
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		int deadline = INT_MAX;
		for (list<Edge*>::iterator eIter = node->out.begin(); eIter != node->out.end(); eIter ++) {
			Edge* edge = *eIter;
			deadline = min(deadline,edge->separationTime);
		}

		node->set_deadline(deadline);

		if (node->out.empty()) {
			cerr << "There exist some nodes without out-going edges!" <<endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// wcets = [1,8]
	for (vector<Node*>::iterator nIter = digraph->node_vec.begin(); nIter != digraph->node_vec.end(); nIter++) {
		Node* node = *nIter;
		node->wcet = 1+rand()%8;
	}
	
	delete[] nodes; 
	return digraph;
}

Digraph** RandomGenerator::generate_digraph_system(int num, int maxNode, double tUtil, bool deadlineProperty) {
	Digraph** digraphs = new Digraph*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numNode = rand()%maxNode+1;
		//int numNode = rand()%21+20;
		//int numNode = rand()%15+1;
		int nScale = sizeof(DIGRAPH_SCALE_FACTOR2)/sizeof(DIGRAPH_SCALE_FACTOR2[0]);
		//int scaleFactor = DIGRAPH_SCALE_FACTOR2[rand()%nScale];
		int scaleFactor = 1;

		Digraph* digraph; 
		//int index = 0;
		//while (true) {
		digraph = generate_one_digraph(i,DIGRAPH_SCALE,numNode);

		int base = calculate_base();

		// generate random edges with random separation time
		for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
			Edge* edge = *iter;
			int separationTime = calculate_separation_time(base);
			edge->set_separation_time(separationTime*DIGRAPH_SCALE*scaleFactor);
		}

		// generate random nodes with random deadline and wcet
		// set arbitrary deadlines
		for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
			Node* node = *iter;

			int minSeparationTime = INT_MAX;
			int maxSeparationTime = INT_MIN;
			for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) {
				minSeparationTime = min(minSeparationTime, (*e)->separationTime);
				maxSeparationTime = max(maxSeparationTime, (*e)->separationTime);
			}

			//int deadline = minSeparationTime + 1.0*(rand()%1000)/1000*(maxSeparationTime-minSeparationTime);
			//int deadline = minSeparationTime + rand()%minSeparationTime + 1;
			int deadline = minSeparationTime + rand()%(2*maxSeparationTime) + 1;
			node->set_deadline(deadline);

			//int wcet = rand()%(deadline/2)+1;
			int wcet = 1.0*(rand()%1000)/1000*minSeparationTime;
			if (wcet == 0) wcet = 1;
			node->set_wcet(wcet);
		}
		
		/// scale wcets to the expected utilization
		// step 1: calculate the linear factor
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		/*
			cout << ++index << endl;

			if (digraph->strongly_connected) {
				break;
			}
			delete digraph;
		}
		*/
		
		digraph->calculate_period_gcd();
		digraph->calculate_linear_factor();
		//cout<<"linear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
		// step 2: rescale wcets
		double factor = util[i]/digraph->linear_factor;
		digraph->scale_wcet(factor);

		/*
		int run = 0;
		while(abs(digraph->linear_factor-util[i])>0.001) {
			double factor = util[i]/digraph->linear_factor;
			digraph->scale_wcet(factor);

			
			// release sccs
			for (vector<Digraph*>::iterator iter = digraph->sccs.begin(); iter != digraph->sccs.end(); iter++) {
				(*iter)->node_vec.clear();
				(*iter)->edge_vec.clear();
				delete *iter;
				*iter = NULL;
			}
			digraph->sccs.clear();
			

			//digraph->generate_strongly_connected_components();
			digraph->calculate_linear_factor();
			cout<<"Run "<<++run<<":\tlinear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
		}
		*/

		digraphs[i] = digraph;
	}

	// reorder digraphs according to the utilization (non-decreasing order)
	// It is possible to considered more.
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];
		order[i] = 0;
		for (int j=0; j<num; j++) {
			Digraph* hdg = digraphs[j];
			//if (hdg->linear_factor < digraph->linear_factor) order[i]++;
			if (hdg->pGCD < digraph->pGCD) order[i]++;
			else if (hdg->pGCD == digraph->pGCD && hdg->linear_factor < digraph->linear_factor) order[i]++;
		}
	}

	Digraph** copy = new Digraph*[num];
	for (int i=0; i<num; i++) copy[i] = digraphs[i];
	for (int i=0; i<num; i++) digraphs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return digraphs;
}

Digraph** RandomGenerator::generate_digraph_system(int num, int minNode, int maxNode, double tUtil, bool deadlineProperty) {
	Digraph** digraphs = new Digraph*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numNode = minNode + rand()%(maxNode-minNode)+1;
		//int numNode = rand()%21+20;
		//int numNode = rand()%15+1;
		int nScale = sizeof(DIGRAPH_SCALE_FACTOR2)/sizeof(DIGRAPH_SCALE_FACTOR2[0]);
		//int scaleFactor = DIGRAPH_SCALE_FACTOR2[rand()%nScale];
		int scaleFactor = 1;

		Digraph* digraph; 
		int index = 0;
		while (true) {
			digraph = generate_one_digraph(i,DIGRAPH_SCALE,numNode);

			int base = calculate_base();

			// generate random edges with random separation time
			for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
				Edge* edge = *iter;
				int separationTime = calculate_separation_time(base);
				edge->set_separation_time(separationTime*DIGRAPH_SCALE*scaleFactor);
			}

			// generate random nodes with random deadline and wcet
			// set arbitrary deadlines
			for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
				Node* node = *iter;

				int minSeparationTime = INT_MAX;
				int maxSeparationTime = INT_MIN;
				for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) {
					minSeparationTime = min(minSeparationTime, (*e)->separationTime);
					maxSeparationTime = max(maxSeparationTime, (*e)->separationTime);
				}

				//int deadline = minSeparationTime + 1.0*(rand()%1000)/1000*(maxSeparationTime-minSeparationTime);
				//int deadline = minSeparationTime + rand()%minSeparationTime + 1;
				//int deadline = minSeparationTime + rand()%(2*maxSeparationTime) + 1;
				//int deadline = minSeparationTime;
				int deadline = minSeparationTime + rand()%(2*maxSeparationTime+1);
				node->set_deadline(deadline);

				//int wcet = rand()%(deadline/2)+1;
				int wcet = 1.0*(rand()%1000)/1000*minSeparationTime;
				if (wcet == 0) wcet = 1;
				node->set_wcet(wcet);
			}

			/// scale wcets to the expected utilization
			// step 1: calculate the linear factor
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			cout << ++index << endl;

			if (!digraph->strongly_connected) {
				delete digraph;
				continue;
			}
			
		
			digraph->calculate_period_gcd();
			digraph->calculate_linear_factor();

			// step 2: rescale wcets
			double factor = util[i]/digraph->linear_factor;
			digraph->scale_wcet(factor);

			break;
		}

		digraphs[i] = digraph;
	}

	// reorder digraphs according to the utilization (non-decreasing order)
	// It is possible to considered more.
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];
		order[i] = 0;
		for (int j=0; j<num; j++) {
			Digraph* hdg = digraphs[j];
			//if (hdg->linear_factor < digraph->linear_factor) order[i]++;
			if (hdg->pGCD < digraph->pGCD) order[i]++;
			else if (hdg->pGCD == digraph->pGCD && hdg->linear_factor < digraph->linear_factor) order[i]++;
		}
	}

	Digraph** copy = new Digraph*[num];
	for (int i=0; i<num; i++) copy[i] = digraphs[i];
	for (int i=0; i<num; i++) digraphs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return digraphs;
}

Digraph** RandomGenerator::generate_digraph_system_without_pointer(int num, int minNode, int maxNode, double tUtil, bool deadlineProperty) {
	Digraph** digraphs = new Digraph*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numNode = minNode + rand()%(maxNode-minNode)+1;
		//int numNode = rand()%21+20;
		//int numNode = rand()%15+1;
		int nScale = sizeof(DIGRAPH_SCALE_FACTOR2)/sizeof(DIGRAPH_SCALE_FACTOR2[0]);
		//int scaleFactor = DIGRAPH_SCALE_FACTOR2[rand()%nScale];
		int scaleFactor = 1;

		Digraph* pDigraph; 
		int index = 0;
		while (true) {
			Digraph digraph = generate_one_digraph_without_pointer(i,DIGRAPH_SCALE,numNode);
			pDigraph = &digraph;
			int base = calculate_base();

			// generate random edges with random separation time
			for (vector<Edge*>::iterator iter = pDigraph->edge_vec.begin(); iter != pDigraph->edge_vec.end(); iter++) {
				Edge* edge = *iter;
				int separationTime = calculate_separation_time(base);
				edge->set_separation_time(separationTime*DIGRAPH_SCALE*scaleFactor);
			}

			// generate random nodes with random deadline and wcet
			// set arbitrary deadlines
			for (vector<Node*>::iterator iter = pDigraph->node_vec.begin(); iter != pDigraph->node_vec.end(); iter++) {
				Node* node = *iter;

				int minSeparationTime = INT_MAX;
				int maxSeparationTime = INT_MIN;
				for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) {
					minSeparationTime = min(minSeparationTime, (*e)->separationTime);
					maxSeparationTime = max(maxSeparationTime, (*e)->separationTime);
				}

				//int deadline = minSeparationTime + 1.0*(rand()%1000)/1000*(maxSeparationTime-minSeparationTime);
				//int deadline = minSeparationTime + rand()%minSeparationTime + 1;
				//int deadline = minSeparationTime + rand()%(2*maxSeparationTime) + 1;
				//int deadline = minSeparationTime;
				int deadline = minSeparationTime + rand()%(2*maxSeparationTime+1);
				node->set_deadline(deadline);

				//int wcet = rand()%(deadline/2)+1;
				int wcet = 1.0*(rand()%1000)/1000*minSeparationTime;
				if (wcet == 0) wcet = 1;
				node->set_wcet(wcet);
			}

			/// scale wcets to the expected utilization
			// step 1: calculate the linear factor
			pDigraph->generate_strongly_connected_components();
			pDigraph->check_strongly_connected();

			cout << ++index << endl;

			if (!pDigraph->strongly_connected) {
				continue;
			}


			pDigraph->calculate_period_gcd();
			pDigraph->calculate_linear_factor();

			// step 2: rescale wcets
			double factor = util[i]/pDigraph->linear_factor;
			pDigraph->scale_wcet(factor);

			break;
		}

		digraphs[i] = pDigraph;
	}

	// reorder digraphs according to the utilization (non-decreasing order)
	// It is possible to considered more.
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];
		order[i] = 0;
		for (int j=0; j<num; j++) {
			Digraph* hdg = digraphs[j];
			//if (hdg->linear_factor < digraph->linear_factor) order[i]++;
			if (hdg->pGCD < digraph->pGCD) order[i]++;
			else if (hdg->pGCD == digraph->pGCD && hdg->linear_factor < digraph->linear_factor) order[i]++;
		}
	}

	Digraph** copy = new Digraph*[num];
	for (int i=0; i<num; i++) copy[i] = digraphs[i];
	for (int i=0; i<num; i++) digraphs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return digraphs;
}

Digraph** RandomGenerator::generate_digraph_system_sccs(int num, int maxNode, double tUtil, bool deadlineProperty) {
	Digraph** digraphs = new Digraph*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numNode = rand()%maxNode+1;
		//int numNode = rand()%21+20;
		//int numNode = rand()%15+1;
		int nScale = sizeof(DIGRAPH_SCALE_FACTOR2)/sizeof(DIGRAPH_SCALE_FACTOR2[0]);
		//int scaleFactor = DIGRAPH_SCALE_FACTOR2[rand()%nScale];
		int scaleFactor = 1;

		Digraph* digraph; 
		int index = 0;
		while (true) {
			digraph = generate_one_digraph(i,DIGRAPH_SCALE,numNode);

			int base = calculate_base();

			// generate random edges with random separation time
			for (vector<Edge*>::iterator iter = digraph->edge_vec.begin(); iter != digraph->edge_vec.end(); iter++) {
				Edge* edge = *iter;

				int separationTime = calculate_separation_time(base);
				edge->set_separation_time(separationTime*1000*scaleFactor);
			}

			// generate random nodes with random deadline and wcet
			// set constrained deadlines
			for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
				Node* node = *iter;

				int minSeparationTime = INT_MAX;
				for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) 
					minSeparationTime = min(minSeparationTime, (*e)->separationTime);

				int deadline = minSeparationTime;
				node->set_deadline(deadline);

				//int wcet = rand()%(deadline/2)+1;
				int wcet = 1.0*(rand()%1000)/1000*deadline;
				if (wcet == 0) wcet = 1;
				node->set_wcet(wcet);
			}

			// set l-MAD deadlines
			if (deadlineProperty) {
				//for (int times=0; times <100; times++) {
					for (vector<Node*>::iterator iter = digraph->node_vec.begin(); iter != digraph->node_vec.end(); iter++) {
						Node* node = *iter;
						//double prob = (double) rand()/RAND_MAX;

						//if (prob <=0.8) {
						//if (true) {
							int deadline = INT_MAX;
							for (list<Edge*>::iterator e = node->out.begin(); e != node->out.end(); e++) {
								Edge* edge = *e;
								int factor = rand()%2+1;
								deadline = min(deadline, edge->snk_node->deadline + factor * edge->separationTime);
							}

							node->set_deadline(deadline);
						//}
					}
				//}
			}
		
			/// scale wcets to the expected utilization
			// step 1: calculate the linear factor
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();
		
			cout << "SCCS Run-" << ++index << endl;

			if (digraph->strongly_connected) {
				break;
			}
			delete digraph;
		}
		
		
		digraph->calculate_period_gcd();
		digraph->calculate_linear_factor();
		//cout<<"linear factor: expected="<<util[i]<<"\tfact="<<digraph->linear_factor<<endl;
		// step 2: rescale wcets
		double factor = util[i]/digraph->linear_factor;
		digraph->scale_wcet(factor);

		digraphs[i] = digraph;
	}

	// reorder digraphs according to the utilization (non-decreasing order)
	// It is possible to considered more.
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Digraph* digraph = digraphs[i];
		order[i] = 0;
		for (int j=0; j<num; j++) {
			Digraph* hdg = digraphs[j];
			//if (hdg->linear_factor < digraph->linear_factor) order[i]++;
			if (hdg->pGCD < digraph->pGCD) order[i]++;
			else if (hdg->pGCD == digraph->pGCD && hdg->linear_factor < digraph->linear_factor) order[i]++;
		}
	}

	Digraph** copy = new Digraph*[num];
	for (int i=0; i<num; i++) copy[i] = digraphs[i];
	for (int i=0; i<num; i++) digraphs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return digraphs;
}

Digraph** RandomGenerator::generate_digraph_system2(int &num, double tUtil) {
	double util = 0;
	int index = 0;
	vector<Digraph*> vec_digraph;
	while (util < tUtil) {
		Digraph* digraph = RandomGenerator::generate_one_digraph3(index, 1);
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		if (!digraph->strongly_connected) {
			delete digraph;
			continue;
		}

		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();

		digraph->calculate_linear_factor();

		cout << "util=" << util << endl;
		cout << "lfac=" << digraph->linear_factor << endl;

		if (digraph->linear_factor+util<= tUtil + 0.01) {
			vec_digraph.push_back(digraph);
			util += digraph->linear_factor;
			index++;

			if (abs(util-tUtil) <= 0.01) break;
		} 
		else
			delete digraph;
	}

	num = vec_digraph.size();
	Digraph** digraphs = new Digraph*[num];
	index=0;
	for (vector<Digraph*>::iterator iter = vec_digraph.begin(); iter != vec_digraph.end(); iter++) {
		digraphs[index] = *iter;
		index++;
	}
	return digraphs;
}

Digraph** RandomGenerator::generate_digraph_system3(int num, double tUtil) {
	Digraph** digraphs = new Digraph*[num];
	double* util = Utility::uniformly_distributed(num,tUtil);

	for (int i=0; i<num; i++) {
		int nScale = sizeof(DIGRAPH_SCALE_FACTOR2)/sizeof(DIGRAPH_SCALE_FACTOR2[0]);
		int scale = DIGRAPH_SCALE_FACTOR2[rand()%nScale];

		Digraph* digraph;
		while (1) {
			digraph = generate_one_digraph4(i,scale,util[i]);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			if (digraph->strongly_connected)
				break;
			delete digraph;
		}

		digraphs[i] = digraph;
	}

	// reorder digraphs according to the scale
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-1; j++) {
			if (digraphs[j]->scale > digraphs[i]->scale) {
				Digraph* temp = digraphs[j];
				digraphs[j] = digraphs[i];
				digraphs[i] = temp;
			}
		}
	}

	return digraphs;
}

Digraph** RandomGenerator::generate_digraph_system4(int &num, double tUtil) {
	double util = 0;
	int index = 0;
	vector<Digraph*> vec_digraph;
	while (util < tUtil) {
		Digraph* digraph = RandomGenerator::generate_DIGRAPH_III(index);
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		if (!digraph->strongly_connected) {
			delete digraph;
			continue;
		}

		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();

		digraph->calculate_linear_factor();

		cout << "util=" << util << endl;
		cout << "lfac=" << digraph->linear_factor << endl;

		if (digraph->linear_factor+util<= tUtil + 0.01) {
			vec_digraph.push_back(digraph);
			util += digraph->linear_factor;
			index++;

			if (abs(util-tUtil) <= 0.01) break;
		} 
		else
			delete digraph;
	}

	num = vec_digraph.size();
	Digraph** digraphs = new Digraph*[num];
	index=0;
	for (vector<Digraph*>::iterator iter = vec_digraph.begin(); iter != vec_digraph.end(); iter++) {
		digraphs[index] = *iter;
		index++;
	}
	return digraphs;
}

Digraph** RandomGenerator::generate_digraph_system_guan(int &num, double tUtil) {
	double util = 0;
	int index = 0;

	vector<Digraph*> vec_digraph;
	int times = 0;

	while (util < tUtil) {
		Digraph* digraph = RandomGenerator::generate_DIGRAPH_VI(index);
		digraph->generate_strongly_connected_components();
		digraph->check_strongly_connected();
		if (!digraph->strongly_connected) {
			delete digraph;
			continue;
		}

		digraph->calculate_period_gcd();
		digraph->calculate_all_gcd();

		digraph->calculate_linear_factor();

		cout << "util=" << util << endl;
		cout << "lfac=" << digraph->linear_factor << endl;

		if (digraph->linear_factor+util<= tUtil + 0.01) {
			vec_digraph.push_back(digraph);
			util += digraph->linear_factor;
			index++;

			if (abs(util-tUtil) <= 0.01) break;
		} 
		else {
			delete digraph;
			times++;
		}

		if (times==5) {
			for (vector<Digraph*>::iterator iter = vec_digraph.begin(); iter != vec_digraph.end(); iter++)
				delete *iter;
			vec_digraph.clear();

			times = 0;
			index = 0;
			util = 0;
		}
	}

	num = vec_digraph.size();
	Digraph** digraphs = new Digraph*[num];
	index=0;
	for (vector<Digraph*>::iterator iter = vec_digraph.begin(); iter != vec_digraph.end(); iter++) {
		digraphs[index] = *iter;
		index++;
	}
	return digraphs;
}

Digraph** RandomGenerator::generate_mixed_digraph_system(int num, double tUtil) {
	int numDigraphI = num/2;
	int numDigraphII = num-numDigraphI;
	
	Digraph** digraphs = new Digraph*[num];

	int size = num;
	// first step: generate DIGRAPH I tasks
	for (int i=0; i<numDigraphI; i++) {
		Digraph* digraph;
		double util;
		if (tUtil < 0) tUtil = 0.001;
		if (size!=1) 
			util = Utility::max_uniformly_distributed(size,tUtil);
		else
			util = tUtil;

		// 100% strongly connected
		int times = 0;
		while(1) {
			digraph = RandomGenerator::generate_DIGRAPH_I(i,util);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			if (digraph->strongly_connected) break;
			times++;
			cout << "times=" << times << endl;

			delete digraph;
		}

		digraph->calculate_linear_factor();
		digraphs[i] = digraph;

		cout << i << "\tExpectedUtil=" << util << "\tRealUtil=" << digraph->linear_factor << endl;
		tUtil -= digraph->linear_factor;
		size--;
	}

	// second step: generate DIGRAPH I tasks
	for (int i=numDigraphI; i<num; i++) {
		Digraph* digraph;
		double util;
		if (tUtil < 0) tUtil = 0.001;
		if (size!=1) 
			util = Utility::max_uniformly_distributed(size,tUtil);
		else
			util = tUtil;
		// 100% strongly connected
		int times = 0;
		while(1) {
			digraph = RandomGenerator::generate_DIGRAPH_II(i,util);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			if (digraph->strongly_connected) break;
			times++;
			cout << "times=" << times << endl;

			delete digraph;
		}

		digraph->calculate_linear_factor();
		digraphs[i] = digraph;

		cout << i << "\tExpectedUtil=" << util << "\tRealUtil=" << digraph->linear_factor << endl;
		tUtil -= digraph->linear_factor;

		size--;
	}

	// reorder digraphs according to the scale
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-1; j++) {
			if (digraphs[j]->scale > digraphs[i]->scale) {
				Digraph* temp = digraphs[j];
				digraphs[j] = digraphs[i];
				digraphs[i] = temp;
			}
		}
	}

	return digraphs;
}


int RandomGenerator::calculate_base() {
	int base = 1;
	int num_ext = rand()%3+1;
	for (int j=0; j<num_ext; j++) base *= DIGRAPH_PERIOD[rand()%3][rand()%2];

	return base;
}

int RandomGenerator::calculate_separation_time(int base) {
	int factor = DIGRAPH_SCALE_FACTOR[rand()%5];

	return factor*base;
}

Stateflow* RandomGenerator::generate_one_stateflow(int index, int scale, int nState, int nTran) {
	Stateflow* sf = new Stateflow(index, scale);

	State** states = new State*[nState];

	for (int i=0; i<nState; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		State* s = new State(name, i);
		states[i] = s;
		sf->add_state(s);
	}

	vector<State*> connected;
	connected.push_back(states[rand()%nState]);
	for (int i=0; i<nState; i++) {
		State* s = states[i];
		
		bool isContained = false;
		for (vector<State*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == s) {
				isContained = true;
				break;
			}
		}

		if (isContained) continue;

		State* s2 = connected.at(rand()%connected.size());
		Transition* t = NULL;

		int b = rand()%2;

		if (b)
			t = new Transition(s,s2);
		else 
			t = new Transition(s2,s);

		sf->add_transition(t);
		connected.push_back(s);
	}

	for (int i=nState-1; i<nTran; i++) {
		State* src = states[rand()%nState];
		State* snk = states[rand()%nState];

		Transition* t = new Transition(src,snk);
		sf->add_transition(t);
	}
	//sf->set_state_number();
	//sf->write_graphviz(cout);
	//int count = 0;
	while (!GraphAlgorithms::isCyclic(sf)) {
		for (int i=0; i<nTran; i++) {
			int reverse = rand()%2;

			if(reverse) {
				Transition* t = sf->trans.at(i);
				
				// reverse the direction of the transition
				t->src->in.push_back(t);
				t->src->out.remove(t);

				t->snk->out.push_back(t);
				t->snk->in.remove(t);

				State* temp = t->src;
				t->src = t->snk;
				t->snk = temp;
			}
		}
		//count ++;
	}
	//cout<<count<<endl;
	//sf->write_graphviz(cout);
	connected.clear();
	delete[] states;
	return sf;
}

Stateflow* RandomGenerator::generate_one_stateflow2(int index, int scale, int nState, int nTran) {
	Stateflow* sf = new Stateflow(index, scale);

	State** states = new State*[nState];

	for (int i=0; i<nState; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		State* s = new State(name, i);
		states[i] = s;
		sf->add_state(s);
	}

	vector<State*> connected;
	connected.push_back(states[rand()%nState]);
	for (int i=0; i<nState; i++) {
		State* s = states[i];
		
		bool isContained = false;
		for (vector<State*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == s) {
				isContained = true;
				break;
			}
		}

		if (isContained) continue;

		State* s2 = connected.at(rand()%connected.size());
		Transition* t = NULL;

		int b = rand()%2;

		if (b)
			t = new Transition(s,s2);
		else 
			t = new Transition(s2,s);

		sf->add_transition(t);
		connected.push_back(s);
	}

	for (int i=nState-1; i<nTran; i++) {
		State* src = states[rand()%nState];
		int count = 0;
		bool exceed = false;
		while (src->out.size() >= 2) {
			src = states[rand()%nState];
			count ++;
			if (count > 100) {
				exceed = true;
				break;
			}
		}
		if (exceed) {nTran = i+1; break;}

		State* snk = states[rand()%nState];

		Transition* t = new Transition(src,snk);
		sf->add_transition(t);
	}
	//sf->set_state_number();
	//sf->write_graphviz(cout);
	//int count = 0;
	while (!GraphAlgorithms::isCyclic(sf)) {
		for (int i=0; i<nTran; i++) {
			int reverse = rand()%2;

			if(reverse) {
				Transition* t = sf->trans.at(i);
				if (t->snk->out.size() >= 1) continue;
				
				// reverse the direction of the transition
				t->src->in.push_back(t);
				t->src->out.remove(t);

				t->snk->out.push_back(t);
				t->snk->in.remove(t);

				State* temp = t->src;
				t->src = t->snk;
				t->snk = temp;
			}
		}
		//count ++;
	}
	//cout<<count<<endl;
	//sf->write_graphviz(cout);
	connected.clear();
	delete[] states;
	return sf;
}

Stateflow* RandomGenerator::generate_one_stateflow3(int index, int scale, int nState) {
	Stateflow* sf = new Stateflow(index, scale);

	State** states = new State*[nState];

	for (int i=0; i<nState; i++) {
		string name = "v"+ Utility::int_to_string(index)+"_"+Utility::int_to_string(i);
		State* s = new State(name, i);
		states[i] = s;
		sf->add_state(s);
	}

	vector<State*> connected;
	for (int i=0; i<nState; i++) {
		State* s = states[i];
		
		bool isContained = false;
		for (vector<State*>::iterator iter=connected.begin(); iter != connected.end(); iter++) {
			if (*iter == s) {
				isContained = true;
				break;
			}
		}

		if (!isContained) connected.push_back(s);

		// 40% with one edge
		// 40% with two edges
		// 10% with three edges
		// 10% with four edges
		double prob = (double) rand()/RAND_MAX;
		int nTran = 0;
		if(prob <= 0.4)      nTran = 1;
		else if(prob <= 0.8) nTran = 2;
		else if(prob <= 0.9) nTran = 3;
		else if(prob <= 1.0) nTran = 4;

		int start = i+1;
		if (connected.size() > i+1) start = i+1+rand()%(connected.size()-i-1);

		nTran = min(nTran, nState-start);
		for (int j=0; j<nTran; j++) {
			State* newState = states[start+j];
			Transition* tran = NULL;

			if (i==0 || (double)rand()/RAND_MAX <= 0.9) 
				tran = new Transition(s,newState);
			else
				tran = new Transition(newState,s);

			sf->add_transition(tran);
		}
	}

	State* src = states[0];
	for (int i=0; i<nState; i++) {
		State* state = states[i];

		// adding backedge from leaf state to src
		if (state->out.size() == 0)  {
			Transition* tran = new Transition(state,src);
			sf->add_transition(tran);
		}
	}

	//sf->set_state_number();
	//sf->write_graphviz(cout);
	//int count = 0;
	while (!GraphAlgorithms::isCyclic(sf)) {
		Transition* t = sf->trans.at(rand()%sf->trans.size());
		if (t->snk->out.size() >= 1) continue;
				
		// reverse the direction of the transition
		t->src->in.push_back(t);
		t->src->out.remove(t);

		t->snk->out.push_back(t);
		t->snk->in.remove(t);

		State* temp = t->src;
		t->src = t->snk;
		t->snk = temp;
	}
	//cout<<count<<endl;
	//sf->write_graphviz(cout);
	connected.clear();
	delete[] states;
	return sf;
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_0(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD3)/sizeof(STATEFLOW_PERIOD3[0]);
	int period_base = STATEFLOW_PERIOD3[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR4)/sizeof(STATEFLOW_FACTOR4[0]);
		int period_factor = STATEFLOW_FACTOR4[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_1(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD4)/sizeof(STATEFLOW_PERIOD4[0]);
	int period_base = STATEFLOW_PERIOD4[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_2(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD4)/sizeof(STATEFLOW_PERIOD4[0]);
	int period_base = STATEFLOW_PERIOD4[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR4)/sizeof(STATEFLOW_FACTOR4[0]);
		int period_factor = STATEFLOW_FACTOR4[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_3(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD3)/sizeof(STATEFLOW_PERIOD3[0]);
	int period_base = STATEFLOW_PERIOD3[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_4(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD3)/sizeof(STATEFLOW_PERIOD3[0]);
	//int period_base = STATEFLOW_PERIOD7[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		/*
		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];
		*/

		t->period = STATEFLOW_PERIOD3[rand()%nPeriod]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_5(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD7)/sizeof(STATEFLOW_PERIOD7[0]);
	//int period_base = STATEFLOW_PERIOD7[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		/*
		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];
		*/

		t->period = STATEFLOW_PERIOD7[rand()%nPeriod]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_6(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD10)/sizeof(STATEFLOW_PERIOD10[0]);

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;
		t->period = STATEFLOW_PERIOD10[rand()%nPeriod]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_7(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD11)/sizeof(STATEFLOW_PERIOD11[0]);

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;
		t->period = STATEFLOW_PERIOD11[rand()%nPeriod]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_8(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;
		t->period = STATEFLOW_PERIOD[rand()%nPeriod]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_exact_analysis_9(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD16)/sizeof(STATEFLOW_PERIOD16[0]);

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;
		t->period = STATEFLOW_PERIOD16[rand()%nPeriod]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_0(Stateflow* sf, int scale) {
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int length = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
		t->period = STATEFLOW_PERIOD[rand()%length]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_1(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
	int period_base = STATEFLOW_PERIOD[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR2)/sizeof(STATEFLOW_FACTOR2[0]);
		int period_factor = STATEFLOW_FACTOR2[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_2(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
	int period_base = STATEFLOW_PERIOD[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_FACTOR3)/sizeof(STATEFLOW_FACTOR3[0]);
		int period_factor = STATEFLOW_FACTOR3[rand()%nFactor];

		t->period = period_base*period_factor*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_3(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int period0 = STATEFLOW_PERIOD5[rand()%nPeriod];
		int period1 = STATEFLOW_PERIOD5[rand()%nPeriod];

		t->period = period0 * period1 * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_4(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
	int period_base = STATEFLOW_PERIOD5[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int period_factor = STATEFLOW_PERIOD5[rand()%nPeriod];

		t->period = period_base * period_factor * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_5(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD6)/sizeof(STATEFLOW_PERIOD6[0]);
	int period_base = STATEFLOW_PERIOD6[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
		int period_factor = STATEFLOW_PERIOD5[rand()%nFactor];

		t->period = period_base * period_factor * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_6(Stateflow* sf, int scale) {
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int length = sizeof(STATEFLOW_PERIOD3)/sizeof(STATEFLOW_PERIOD3[0]);
		t->period = STATEFLOW_PERIOD3[rand()%length]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_7(Stateflow* sf, int scale) {
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int length = sizeof(STATEFLOW_PERIOD4)/sizeof(STATEFLOW_PERIOD4[0]);
		t->period = STATEFLOW_PERIOD4[rand()%length]*scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_8(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD8)/sizeof(STATEFLOW_PERIOD8[0]);
	int period_base = STATEFLOW_PERIOD8[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
		int period_factor = STATEFLOW_PERIOD5[rand()%nFactor];

		t->period = period_base * period_factor * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_9(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD9)/sizeof(STATEFLOW_PERIOD9[0]);
	int period_base = STATEFLOW_PERIOD9[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
		int period_factor = STATEFLOW_PERIOD5[rand()%nFactor];

		t->period = period_base * period_factor * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_10(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD12)/sizeof(STATEFLOW_PERIOD12[0]);
	int period_base = STATEFLOW_PERIOD12[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
		int period_factor = STATEFLOW_PERIOD5[rand()%nFactor];

		t->period = period_base * period_factor * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_11(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD13)/sizeof(STATEFLOW_PERIOD13[0]);
	int period_base = STATEFLOW_PERIOD13[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
		int period_factor = STATEFLOW_PERIOD5[rand()%nFactor];

		t->period = period_base * period_factor * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_12(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD14)/sizeof(STATEFLOW_PERIOD14[0]);
	int period_base = STATEFLOW_PERIOD14[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
		int period_factor = STATEFLOW_PERIOD5[rand()%nFactor];

		t->period = period_base * period_factor * scale;

		t->priority = k++;
	}
}

void RandomGenerator::setup_periods_for_one_stateflow_for_approximate_analysis_13(Stateflow* sf, int scale) {
	int nPeriod = sizeof(STATEFLOW_PERIOD17)/sizeof(STATEFLOW_PERIOD17[0]);
	int period_base = STATEFLOW_PERIOD17[rand()%nPeriod];

	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int nFactor = sizeof(STATEFLOW_PERIOD5)/sizeof(STATEFLOW_PERIOD5[0]);
		int period_factor = STATEFLOW_PERIOD5[rand()%nFactor];

		t->period = period_base * period_factor * scale;

		//if (t->period <=0 ) t->period = 2000000000;
		//	cout << "heeee" <<endl;

		t->priority = k++;
	}
}

void RandomGenerator::setup_wcets_for_one_stateflow(Stateflow* sf, double util) {
	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int wcet = (int) (util*t->period);
		if (wcet == 0) wcet = 1;

		//int wcet = rand()%5+1;
		t->wcet = wcet;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	cout << "n_time_instance=" << sf->n_time_instance << endl;
	cout << "hyperperiod=" << sf->hyperperiod << endl;
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	// step 2: rescale wcets
	int run = 0;
	while (abs(sf->lfac-util)>0.0001) {
		double factor = util/sf->lfac;
		if (sf->lfac == 0) 
			cout << "here0" <<endl;
		sf->scale_wcet(factor);
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;
		if (sf->lfac == 0) 
			cout << "here1" <<endl;
	}
}

void RandomGenerator::setup_wcets_for_one_stateflow(Stateflow* sf, double util, double& diffUtil) {
	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int wcet = (int) (util*t->period);
		if (wcet == 0) wcet = 1;

		//int wcet = rand()%5+1;
		t->wcet = wcet;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	cout << "n_time_instance=" << sf->n_time_instance << endl;
	cout << "hyperperiod=" << sf->hyperperiod << endl;
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	// step 2: rescale wcets for only once
	double factor = util/sf->lfac;
	sf->scale_wcet(factor);
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	diffUtil = util-sf->lfac;

	cout<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac << "\tdiffUtil=" << diffUtil <<endl;
}

void RandomGenerator::setup_wcets_for_one_stateflow2(Stateflow* sf, double util) {
	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int wcet = (int) (util*t->period);
		if (wcet == 0) wcet = 1;

		//int wcet = rand()%5+1;
		t->wcet = wcet;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	cout << "n_time_instance=" << sf->n_time_instance << endl;
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	cout<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	/*
	// step 2: rescale wcets
	int run = 0;
	while (abs(sf->lfac-util)>0.0001) {
		double factor = util/sf->lfac;
		sf->scale_wcet(factor);
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;
	}
	*/
}

Stateflow* RandomGenerator::generate_one_stateflow_with_util(int index, int scale, double util, int nState, int nTran) {
	Stateflow* sf = generate_one_stateflow(index,scale,nState, nTran);

	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;

		int length = sizeof(STATEFLOW_PERIOD)/sizeof(STATEFLOW_PERIOD[0]);
		t->period = STATEFLOW_PERIOD[rand()%length]*scale;

		int wcet = (int) (t->period*(rand()%100/100.0));
		if (wcet == 0) wcet = rand()%5+1;

		//int wcet = rand()%5+1;
		t->wcet = wcet;

		t->priority = k++;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	// step 2: rescale wcets
	int run = 0;
	while (abs(sf->lfac-util)>0.001) {
		double factor = util/sf->lfac;
		sf->scale_wcet(factor);
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		//cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;
	}

	return sf;
}

Stateflow* RandomGenerator::generate_one_stateflow_with_util2(int index, int scale, double util, int nState, int nTran) {
	Stateflow* sf = generate_one_stateflow(index,scale,nState, nTran);

	int base_length = sizeof(STATEFLOW_BASE)/sizeof(STATEFLOW_BASE[0]);
	int base = STATEFLOW_BASE[rand()%base_length];
	//double* tran_util = Utility::uniformly_distributed(numTran,1);
	int k=0;
	for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
		Transition* t = *iter;
		int factor_length = sizeof(STATEFLOW_FACTOR)/sizeof(STATEFLOW_FACTOR[0]);
		int factor = STATEFLOW_FACTOR[rand()%factor_length];
		t->period = base*factor;

		//int wcet = (int) (t->period*(rand()%100/100.0));
		//if (wcet == 0) wcet = rand()%5+1;
		int wcet = factor*(rand()%5+1);
		t->wcet = wcet;

		t->priority = k++;
	}

	/// scale wcets to the expected utilization
	// step 1: calculate the linear factor
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	//sf->generate_precise_digraph();
	//cout<<sf->lfac<<"\t"<<sf->precise_digraph->linear_factor<<endl;

	// step 2: rescale wcets
	int run = 0;
	while (abs(sf->lfac-util)>0.001) {
		double factor = util/sf->lfac;
		sf->scale_wcet(factor);
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		//cout<<"Run:"<<run++<<"\tlinear factor: expected="<<util<<"\tfact="<<sf->lfac<<endl;
	}

	return sf;
}

Stateflow** RandomGenerator::generate_stateflow_system(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	for (int i=0; i<num; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system2(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	for (int i=0; i<num; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow_with_util2(i,STATEFLOW_SCALE, util[i], numState,numTran);;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system3(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	double* util = Utility::uniformly_distributed(num, tUtil);

	int bound = 2*num/3;

	for (int i=0; i<bound; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow_with_util(i,STATEFLOW_SCALE, util[i], numState,numTran);
	}

	for (int i=bound; i<num; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow_with_util2(i,STATEFLOW_SCALE, util[i], numState,numTran);
	}

	delete[] util;
	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system4(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	for (int i=0; i<num-1; i++) {
		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num-1];
	for (int i=0; i<num-1; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num-1; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num-1];
	for (int i=0; i<num-1; i++) copy[i] = sfs[i];
	for (int i=0; i<num-1; i++) sfs[order[i]] = copy[i];

	double* util = Utility::uniformly_distributed(num-1, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num-1; i++) {
		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	string name = "Input\\OneSpecial.dot";
	const char *p = name.c_str();

	Stateflow* sf = FileReader::ReadOneStateflow(1,p);
	sf->calculate_gcd();
	sf->calculate_hyperperiod();
	sf->set_state_number();
	sf->generate_rbf_time_instances();
	sf->generate_rbfs();
	sf->generate_exec_req_matrix();
	//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
	sf->calculate_linear_factor();

	sfs[num-1] = sf;

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system5(int num, int maxState, double tUtil) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = rand()%(maxState-1)+2;
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system6(int num, int maxState, double tUtil, int period_choice) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.45) numState = 3;
		else if (r<0.55) numState = 4;
		else if (r<0.65) numState = 5;
		else if (r<0.7) numState = 6;
		else if (r<0.75) numState = 7;
		else if (r<0.8) numState = 8;
		else if (r<0.85) numState = 9;
		else if (r<0.9) numState = 10;
		else if (r<0.94) numState = 11;
		else if (r<0.97) numState = 12;
		else if (r<0.99) numState = 13;
		else if (r<=1) numState = 14;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		cout<<"numState="<<numState<<endl;
		
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system7(int num, int maxState, double tUtil, int period_choice) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.45) numState = 3;
		else if (r<0.55) numState = 4;
		else if (r<0.65) numState = 5;
		else if (r<0.7) numState = 6;
		else if (r<0.75) numState = 7;
		else if (r<0.8) numState = 8;
		else if (r<0.85) numState = 9;
		else if (r<0.9) numState = 10;
		else if (r<0.94) numState = 11;
		else if (r<0.97) numState = 12;
		else if (r<0.99) numState = 13;
		else if (r<=1) numState = 14;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		cout<<"numState="<<numState<<endl;
		
		int numTran = calculate_num_edge(numState);

		sfs[i] = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system8(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.45) numState = 3;
		else if (r<0.55) numState = 4;
		else if (r<0.65) numState = 5;
		else if (r<0.7) numState = 6;
		else if (r<0.75) numState = 7;
		else if (r<0.8) numState = 8;
		else if (r<0.85) numState = 9;
		else if (r<0.9) numState = 10;
		else if (r<0.94) numState = 11;
		else if (r<0.97) numState = 12;
		else if (r<0.99) numState = 13;
		else if (r<=1) numState = 14;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		cout<<"numState="<<numState<<endl;
		
		int numTran = calculate_num_edge(numState);

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			/*
			if (digraph->strongly_connected && r<scc_probability) {
				delete digraph;
				break;
			}
			else if (!digraph->strongly_connected && r<1-scc_probability) {
				delete digraph;
				break;
			}
			else {
				delete digraph;
				delete sf;
			}
			*/
			delete digraph;
			delete sf;

		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);
	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/
	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_exact_analysis(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;
		// number of states = {20%--1, 15%--2, 
		// 15%--3, 10%--4, 10%--5, 10%--6, 5%--7, 
		// 5%--8, 5%--9, 5%--10}, which makes the smaller number 
		// of states has higher proportion
		/*
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.5) numState = 3;
		else if (r<0.6) numState = 4;
		else if (r<0.7) numState = 5;
		else if (r<0.8) numState = 6;
		else if (r<0.85) numState = 7;
		else if (r<0.9) numState = 8;
		else if (r<0.95) numState = 9;
		else if (r<=1) numState = 10;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		*/

		
		// number of states = {30%--1, 20%--2, 20%--3, 20%--4, 10%--5}
		if (r<0.3) numState = 1;
		else if (r<0.5) numState = 2;
		else if (r<0.7) numState = 3;
		else if (r<0.9) numState = 4;
		else if (r<=1) numState = 5;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		
		
		int numTran = calculate_num_edge(numState);
		//numTran = (numTran+numState)/2;
		cout<<"numState="<<numState<<"\t"<<"numTran"<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;

		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis4)
			setup_periods_for_one_stateflow_for_exact_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;

		/*
		// fist configure
		// number of states = {20%--1, 15%--2, 
		// 15%--3, 10%--4, 10%--5, 10%--6, 5%--7, 
		// 5%--8, 5%--9, 5%--10}, which makes the smaller number 
		// of states has higher proportion
		if (r<0.2) numState = 1;
		else if (r<0.35) numState = 2;
		else if (r<0.5) numState = 3;
		else if (r<0.6) numState = 4;
		else if (r<0.7) numState = 5;
		else if (r<0.8) numState = 6;
		else if (r<0.85) numState = 7;
		else if (r<0.9) numState = 8;
		else if (r<0.95) numState = 9;
		else if (r<=1) numState = 10;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		*/

		
		// second configure
		// number of states = {15%--1, 15%--2, 
		// 10%--3, 10%--4, 10%--5, 5%---6, 5%---7, 
		// 5%---8, 5%---9, 5%---10, 5%---11, 4%---12, 
		// 3%---13, 2%---14, 1%---15}

		if (r<0.15) numState = 1;
		else if (r<0.30) numState = 2;
		else if (r<0.40) numState = 3;
		else if (r<0.50) numState = 4;
		else if (r<0.60) numState = 5;
		else if (r<0.65) numState = 6;
		else if (r<0.70) numState = 7;
		else if (r<0.75) numState = 8;
		else if (r<0.80) numState = 9;
		else if (r<0.85) numState = 10;
		else if (r<0.9) numState = 11;
		else if (r<0.94) numState = 12;
		else if (r<0.97) numState = 13;
		else if (r<=0.99) numState = 14;
		else if (r<=1) numState = 15;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		
		
		// third configure
		//numState = rand()%maxState+1;

		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran"<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;
		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis2(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;

		// state number configure
		// 20% with 1 state
		// 20% with 2 states 
		// 20% with 3 states
		// 20% with 5 states
		// 10% with 10 states
		// 5% with 15 states
		// 5% with 20 states

		if (r<0.20) numState = 1;
		else if (r<0.40) numState = 2;
		else if (r<0.60) numState = 3;
		else if (r<0.80) numState = 5;
		else if (r<0.90) numState = 10;
		else if (r<0.95) numState = 15;
		else if (r<=1) numState = 20;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}

		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran="<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;
		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis3)
			setup_periods_for_one_stateflow_for_approximate_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis4)
			setup_periods_for_one_stateflow_for_approximate_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis5)
			setup_periods_for_one_stateflow_for_approximate_analysis_5(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis3(int num, int numState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	for (int i=0; i<num; i++) {
		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran="<<numTran<<endl;

		Stateflow* sf;

		sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis3)
			setup_periods_for_one_stateflow_for_approximate_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis4)
			setup_periods_for_one_stateflow_for_approximate_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis5)
			setup_periods_for_one_stateflow_for_approximate_analysis_5(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis6)
			setup_periods_for_one_stateflow_for_approximate_analysis_6(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis7)
			setup_periods_for_one_stateflow_for_approximate_analysis_7(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis8)
			setup_periods_for_one_stateflow_for_approximate_analysis_8(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis9)
			setup_periods_for_one_stateflow_for_approximate_analysis_9(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis13)
			setup_periods_for_one_stateflow_for_approximate_analysis_13(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis4(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;

		// state number configure
		// 20% with 1 state
		// 20% with 2 states 
		// 20% with 3 states
		// 20% with 5 states
		// 10% with 10 states
		// 5% with 15 states
		// 5% with 20 states
		/*
		if (r<0.20) numState = 1;
		else if (r<0.40) numState = 2;
		else if (r<0.60) numState = 3;
		else if (r<0.80) numState = 5;
		else if (r<0.90) numState = 10;
		else if (r<0.95) numState = 15;
		else if (r<=1) numState = 20;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		*/
		// state number configure
		// 20% with 1 state
		// 20% with 5 states 
		// 20% with 10 states
		// 20% with 15 states
		// 10% with 20 states
		// 5% with 25 states
		// 5% with 30 states
		/*
		if (r<0.20) numState = 10;
		else if (r<0.40) numState = 20;
		else if (r<0.60) numState = 25;
		else if (r<0.80) numState = 50;
		else if (r<0.90) numState = 75;
		//else if (r<0.95) numState = 100;
		else if (r<=1) numState = 100;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}
		*/

		if (r<0.20) numState = 5;
		else if (r<0.40) numState = 10;
		else if (r<0.60) numState = 15;
		else if (r<0.80) numState = 20;
		else if (r<0.90) numState = 25;
		//else if (r<0.95) numState = 100;
		else if (r<=1) numState = 50;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}

		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran="<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;
		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis3)
			setup_periods_for_one_stateflow_for_approximate_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis4)
			setup_periods_for_one_stateflow_for_approximate_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis5)
			setup_periods_for_one_stateflow_for_approximate_analysis_5(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis6)
			setup_periods_for_one_stateflow_for_approximate_analysis_6(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis7)
			setup_periods_for_one_stateflow_for_approximate_analysis_7(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis8)
			setup_periods_for_one_stateflow_for_approximate_analysis_8(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis9)
			setup_periods_for_one_stateflow_for_approximate_analysis_9(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis5(int num, int numState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	for (int i=0; i<num; i++) {
		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran="<<numTran<<endl;

		Stateflow* sf;

		sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis3)
			setup_periods_for_one_stateflow_for_approximate_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis4)
			setup_periods_for_one_stateflow_for_approximate_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis5)
			setup_periods_for_one_stateflow_for_approximate_analysis_5(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed2(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis6(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;

		// state number configure
		// 20% with 1 state
		// 20% with 5 states 
		// 20% with 10 states
		// 20% with 15 states
		// 10% with 20 states
		// 5% with 25 states
		// 5% with 30 states

		if (r<0.20) numState = 1;
		else if (r<0.40) numState = 5;
		else if (r<0.60) numState = 10;
		else if (r<0.80) numState = 15;
		else if (r<0.90) numState = 20;
		else if (r<0.95) numState = 25;
		else if (r<=1) numState = 30;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}

		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran="<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;
		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis3)
			setup_periods_for_one_stateflow_for_approximate_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis4)
			setup_periods_for_one_stateflow_for_approximate_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis5)
			setup_periods_for_one_stateflow_for_approximate_analysis_5(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis6)
			setup_periods_for_one_stateflow_for_approximate_analysis_6(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis7)
			setup_periods_for_one_stateflow_for_approximate_analysis_7(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis7(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;

		// state number configure
		// 20% with 1 state
		// 20% with 2 states 
		// 20% with 3 states
		// 20% with 5 states
		// 10% with 10 states
		// 5% with 15 states
		// 5% with 20 states
		
		if (r<0.20) numState = 1;
		else if (r<0.40) numState = 2;
		else if (r<0.60) numState = 3;
		else if (r<0.80) numState = 5;
		else if (r<0.90) numState = 10;
		else if (r<0.95) numState = 15;
		else if (r<=1) numState = 20;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}

		int numTran = calculate_num_edge(numState);
		cout<<"numState="<<numState<<"\t"<<"numTran="<<numTran<<endl;

		Stateflow* sf;

		while (true) {
			sf = generate_one_stateflow2(i,STATEFLOW_SCALE,numState, numTran);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			if (r<1-scc_probability) {
				delete digraph;
				break;
			}

			delete digraph;
			delete sf;
		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis4)
			setup_periods_for_one_stateflow_for_exact_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis5)
			setup_periods_for_one_stateflow_for_exact_analysis_5(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis6)
			setup_periods_for_one_stateflow_for_exact_analysis_6(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis7)
			setup_periods_for_one_stateflow_for_exact_analysis_7(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis8)
			setup_periods_for_one_stateflow_for_exact_analysis_8(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis9)
			setup_periods_for_one_stateflow_for_exact_analysis_9(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis3)
			setup_periods_for_one_stateflow_for_approximate_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis4)
			setup_periods_for_one_stateflow_for_approximate_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis5)
			setup_periods_for_one_stateflow_for_approximate_analysis_5(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis6)
			setup_periods_for_one_stateflow_for_approximate_analysis_6(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis7)
			setup_periods_for_one_stateflow_for_approximate_analysis_7(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis8)
			setup_periods_for_one_stateflow_for_approximate_analysis_8(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis9)
			setup_periods_for_one_stateflow_for_approximate_analysis_9(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}

int RandomGenerator::calculate_num_edge(int numNode) {
	int divisor = 1;
    for(int j=0; j<numNode; j++) divisor = divisor + rand()%2;
    	
    int s = ceil(1.0*(rand()%numNode)/divisor); 
    s = min(s, numNode + 1);
    	
    int numEdge = numNode;
    for(int j=0; j<s; j++) numEdge = numEdge + rand()%(numNode + 1);
		
    return numEdge;
}

Stateflow** RandomGenerator::generate_stateflow_system_for_approximate_analysis_with_periodicity(int num, int maxState, double tUtil, int period_choice, double scc_probability) {
	Stateflow** sfs = new Stateflow*[num];

	int sIndex = -1;
	bool isExisted = false;
	//if (rand()%100<80) isExisted = true;

	if (isExisted) sIndex = rand()%num;

	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		int numState = 0;
		double r = (double)rand()/RAND_MAX;

		if (r<0.20) numState = 10;
		else if (r<0.40) numState = 15;
		else if (r<0.60) numState = 20;
		else if (r<0.80) numState = 25;
		else if (r<0.90) numState = 40;
		//else if (r<0.95) numState = 100;
		else if (r<=1) numState = 50;
		else {
			cerr<<"Error random double "<<r<<endl;
			exit(EXIT_FAILURE);
		}

		//numState = rand()%maxState+1;

		//int numTran = calculate_num_edge(numState);
		//cout<<"numState="<<numState<<"\t"<<"numTran="<<numTran<<endl;

		Stateflow* sf;
		int sccNum = 1;

		while (true) {
			sf = generate_one_stateflow3(i,STATEFLOW_SCALE,numState);
			Digraph* digraph = GraphAlgorithms::generate_simple_digraph(sf);
			digraph->generate_strongly_connected_components();
			digraph->check_strongly_connected();

			//r = (double)rand()/RAND_MAX;

			if (digraph->strongly_connected) {
				delete digraph;
				break;
			}

			/*
			if (r<1-scc_probability) {
				delete digraph;
				break;
			}
			*/

			//cout << "SCC Run: " << sccNum++ << "\tnumState="<<numState<<"\t"<<"numTran="<<numTran<<endl;
			cout << "SCC Run: " << sccNum++ << "\tnumState="<<numState<<endl;
			//numTran = (sccNum/10+1) * numTran;

			delete digraph;
			delete sf;
		}

		sfs[i] = sf;
		// set up the periods for the stateflow
		// in general, period_choice should be ApproximateAnalysis1
		if (period_choice == ExactAnalysis0)
			setup_periods_for_one_stateflow_for_exact_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis1)
			setup_periods_for_one_stateflow_for_exact_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis2)
			setup_periods_for_one_stateflow_for_exact_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ExactAnalysis3)
			setup_periods_for_one_stateflow_for_exact_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis0)
			setup_periods_for_one_stateflow_for_approximate_analysis_0(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis1)
			setup_periods_for_one_stateflow_for_approximate_analysis_1(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis2)
			setup_periods_for_one_stateflow_for_approximate_analysis_2(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis3)
			setup_periods_for_one_stateflow_for_approximate_analysis_3(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis4)
			setup_periods_for_one_stateflow_for_approximate_analysis_4(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis5)
			setup_periods_for_one_stateflow_for_approximate_analysis_5(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis6)
			setup_periods_for_one_stateflow_for_approximate_analysis_6(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis7)
			setup_periods_for_one_stateflow_for_approximate_analysis_7(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis8)
			setup_periods_for_one_stateflow_for_approximate_analysis_8(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis9)
			setup_periods_for_one_stateflow_for_approximate_analysis_9(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis10)
			setup_periods_for_one_stateflow_for_approximate_analysis_10(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis11)
			setup_periods_for_one_stateflow_for_approximate_analysis_11(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis12)
			setup_periods_for_one_stateflow_for_approximate_analysis_12(sfs[i], STATEFLOW_SCALE);
		else if (period_choice == ApproximateAnalysis13)
			setup_periods_for_one_stateflow_for_approximate_analysis_13(sfs[i], STATEFLOW_SCALE);
		else {
			cerr<<"Error Period Choice = "<<period_choice<<endl;
			exit(EXIT_FAILURE);
		}
		
	}

	double* util = Utility::uniformly_distributed(num, tUtil);

	/*
	// reorder utilizations by non-decreasing order
	for (int i=0; i<num; i++) {
		for (int j=0; j<num-i-1; j++) {
			if (util[j]>util[j+1]) {
				double temp = util[j];
				util[j] = util[j+1];
				util[j+1] = temp;
			}
		}
	}
	*/

	// set up the wcets for the stateflow based on the expected utilization
	for (int i=0; i<num; i++) {
		if (i==sIndex) continue;

		setup_wcets_for_one_stateflow(sfs[i], util[i]);
		//setup_wcets_for_one_stateflow(sfs[i], tUtil/num);
	}

	// for the special one
	if(sIndex != -1) {
		cout<<sIndex<<endl;
		string name = "Input\\OneSpecial.dot";
		const char *p = name.c_str();

		Stateflow* sf = FileReader::ReadOneStateflow(1,p);
		// set wcets
		for (vector<Transition*>::iterator iter = sf->trans.begin(); iter != sf->trans.end(); iter++) {
			Transition* t = *iter;
			if (t->period==100000) t->wcet = 10000+rand()%5000;
			if (t->period==20000) {
				int wcet = util[sIndex]*t->period;
				if (wcet == 0) wcet = rand()%5+1;
				t->wcet = wcet;
			}
		}

		sf->calculate_gcd();
		sf->calculate_hyperperiod();
		sf->set_state_number();
		sf->generate_rbf_time_instances();
		sf->generate_rbfs();
		sf->generate_exec_req_matrix();
		//Utility::output_matrix(sf->exec_req_matrix,sf->n_state,sf->n_state);
		sf->calculate_linear_factor();

		sfs[sIndex] = sf;
	}

	// reorder stateflows according to the gcd (an estimate on the smallest deadline)
	int* order = new int[num];
	for (int i=0; i<num; i++) {
		Stateflow* sf = sfs[i];
		order[i] = 0;
		for ( int j=0; j<num; j++) {
			Stateflow* hsf = sfs[j];
			if(hsf->gcd < sf->gcd) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod<sf->hyperperiod) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state < sf->n_state) order[i]++;
			else if (hsf->gcd == sf->gcd && hsf->hyperperiod == sf->hyperperiod && hsf->n_state == sf->n_state && i<j) order[i]++;
		}
	}

	Stateflow** copy = new Stateflow*[num];
	for (int i=0; i<num; i++) copy[i] = sfs[i];
	for (int i=0; i<num; i++) sfs[order[i]] = copy[i];

	delete[] util;
	delete[] order;
	delete[] copy;

	return sfs;
}