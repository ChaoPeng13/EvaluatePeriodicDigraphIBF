#include "Stateflow.h"

#include "MaxPlusAlgebra.h"
#include "GraphAlgorithms.h"

#include "Timer.h"

double Stateflow::tDiff = 0;
double Stateflow::tCalWithoutPeriodicity = 0;
double Stateflow::tCalWithPeriodicity = 0;

Stateflow::~Stateflow() {
	//cout<<"Stateflow Destruction Start"<<endl;
	state_index.clear();
	index_state.clear();
	tran_index.clear();
	index_tran.clear();
	// release states
	for (vector<State*>::iterator iter = states.begin(); iter != states.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	states.clear();
		
	// release transitions
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	trans.clear();

	rbf_time_instances.clear();
	index_time.clear();
	time_index.clear();

	//rbf_hyperperiod.clear();
	//ibf_hyperperiod.clear();
		
	// delete execution request matrix
	for(map<int, double**>::iterator iter = exec_req_matrix_power.begin(); iter != exec_req_matrix_power.end(); iter++) {
		double** temp = iter->second;
		for (int i=0; i<n_state; i++) {
			delete[] temp[i];
		}
		delete[] temp;
	}
	exec_req_matrix_power.clear();

	if (simple_digraph != NULL) delete simple_digraph;
	if (precise_digraph != NULL) delete precise_digraph;
	if (sc_precise_digraph != NULL) delete sc_precise_digraph;
	if (reachability_digraph != NULL) delete reachability_digraph;
	//if (exec_digraph != NULL) delete exec_digraph;

	// release rbf_vec
	for (vector<StartFinishTime*>::iterator iter = rbf_vec.begin(); iter != rbf_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	rbf_vec.clear();

	// release ibf_vec
	for (vector<StartFinishTime*>::iterator iter = ibf_vec.begin(); iter != ibf_vec.end(); iter++) {
		delete *iter;
		*iter = NULL;
	}
	ibf_vec.clear();

	// release abstract request function trees
	for (map<int, map<int,AbstractRequestFunctionTree>>::iterator mmIter = mmARFT.begin(); mmIter != mmARFT.end(); mmIter++) {
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

	// the four-dimension matrix rbfs has been released in the ends of generate_rbfs or generate_ibfs
	// release rbfsf
	if (rbfsf != NULL) {
		for (int s=0; s<n_time_instance; s++) {
			delete[] rbfsf[s];
		}
		delete[] rbfsf;
	}

	// release rbfif
	if (rbfif != NULL) {
		for (int i=0; i<n_state; i++) {
			delete[] rbfif[i];
		}
		delete[] rbfif;
	}

	// relase rbfjs
	if (rbfjs != NULL) {
		for (int j=0; j<n_state; j++) {
			delete[] rbfjs[j];
		}
		delete[] rbfjs;
	}

	// release ibfjs
	if(ibfjs != NULL) {
		for (int j=0; j<n_state; j++) {
			delete[] ibfjs[j];
		}
		delete[] ibfjs;
	}

	//cout<<"Staeflow End"<<endl;
}

void Stateflow::add_state(State* state) {
	states.push_back(state);
	state_index[state] = iState;
	index_state[iState] = state;
	iState++;
}

void Stateflow::add_transition(Transition* tran) {
	trans.push_back(tran);
	tran_index[tran] = iTran;
	index_tran[iTran] = tran;

	tran->src->out.push_back(tran);
	tran->snk->in.push_back(tran);

	iTran++;
}

bool Stateflow::less_compare_trans(const Transition* t1, const Transition* t2) {
	return t1->priority < t2->priority;
}

void Stateflow::prepare() {
	// TODO:
}

void Stateflow::calculate_gcd() {
	int _gcd = 0;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		_gcd = Utility::math_gcd(_gcd, tran->period);
	}

	gcd = _gcd;
}

void Stateflow::calculate_t_gcd() {
	int _gcd = 0;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		_gcd = Utility::math_gcd(_gcd, tran->period);
		_gcd = Utility::math_gcd(_gcd, tran->wcet);
	}
	t_gcd = _gcd;
}

void Stateflow::calculate_hyperperiod() {
	int _lcm = 1;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		_lcm = Utility::math_lcm(_lcm, tran->period);
	}
	hyperperiod = _lcm;
}

void Stateflow::calculate_deadlines() {
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		State* snk = (*iter)->snk;
		int deadline = (*iter)->period;
		for (list<Transition*>::iterator iter2 = snk->out.begin(); iter2 != snk->out.end(); iter2++) {
			deadline = min(deadline,(*iter2)->period);
		}
		(*iter)->deadline = deadline;
	}
}

void Stateflow::set_state_number() {
	n_state = states.size();
}

void Stateflow::generate_rbf_time_instances() {
	rbf_time_instances.insert(0);

	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		for (int i=tran->period; i<hyperperiod; i+=tran->period)
			rbf_time_instances.insert(i);
	}

	rbf_time_instances.insert(hyperperiod);
	int index = 0;
	for (set<int>::iterator iter = rbf_time_instances.begin(); iter != rbf_time_instances.end();  iter++) {
		index_time[index] = *iter;
		time_index[*iter] = index;
		index++;
	}

	n_time_instance = rbf_time_instances.size();
}

void Stateflow::generate_ibf_time_instances() {
	// TODO: It needs to think more
}

void Stateflow::generate_rbfs() {
	// Sort the transitions by the priorities
	//sort(trans.begin(),trans.end(),less_compare_trans);
	// allocate memory for rbf_{i,j}[s,f)
	//if (rbfs == NULL) {
		rbfs = new double***[n_state];
		for (int i=0; i<n_state; i++) rbfs[i] = new double**[n_state];
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) rbfs[i][j] = new double*[n_time_instance];
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) for (int k=0; k<n_time_instance; k++) rbfs[i][j][k] = new double[n_time_instance];
	//}

	//cout << n_state*n_state*n_time_instance*n_time_instance << endl;

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++)
		for (int s=0; s<n_time_instance; s++) for (int f=0; f<n_time_instance; f++) {
			if (i==j && s<=f) rbfs[i][j][s][f] = 0;
			else rbfs[i][j][s][f] = NEG_INFINITY;
		}
	
	// Dynamic Programming
	for (int len=1; len<n_time_instance; len++) {
		
		for (int s=0; s<n_time_instance-len; s++) {
			int f = s+len; // index of finish time
			int prev_ft = index_time[f-1]; // previous finish time

			//initialize
			for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
				if (len==1) {
					if (i==j) rbfs[i][j][s][f]=0;
				} else rbfs[i][j][s][f] = rbfs[i][j][s][f-1];
			}

			// calculate rbfs[i][j][s][f] based on the assumption we have known all the rbf[x][y][s][s..f-1].

			for (vector<Transition*>::iterator t_iter = trans.begin(); t_iter != trans.end(); t_iter++) {
				Transition* tran = *t_iter;
				int srcIndx = state_index[tran->src];
				int snkIndx = state_index[tran->snk];

				if(prev_ft%tran->period == 0) { // this transition might happen
					for (vector<State*>::iterator s_iter = states.begin(); s_iter != states.end(); s_iter++) {
						State* state = *s_iter;
						int sIndx = state_index[state];

						if (len == 1 && sIndx == srcIndx)
							rbfs[sIndx][snkIndx][s][f] = max(rbfs[sIndx][snkIndx][s][f], (double)tran->wcet);
						else if (len >1)
							rbfs[sIndx][snkIndx][s][f] = max(rbfs[sIndx][snkIndx][s][f], rbfs[sIndx][srcIndx][s][f-1]+tran->wcet);
					}
				}
			}
		}
	}

	// set rbfsf
	if (rbfsf == NULL) { // allocate the memory for rbf[s,f)
		rbfsf = new double*[n_time_instance];
		for (int s=0; s<n_time_instance; s++) rbfsf[s] = new double[n_time_instance]; 
		for (int s=0; s<n_time_instance; s++) for (int f=0; f<n_time_instance; f++) {
			if ( s <= f) rbfsf[s][f] = 0.0;
			else rbfsf[s][f] = NEG_INFINITY;
		}
	}

	for (int s=0; s<n_time_instance; s++) for (int f=0; f<n_time_instance; f++) {
		double temp = NEG_INFINITY;
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
			temp = max(temp, rbfs[i][j][s][f]);
		}
		rbfsf[s][f] = temp;
	}

	// set rbfif
	if (rbfif == NULL) { // allocate the memory for rbf_{i}[0,f)
		rbfif = new double*[n_state];
		for (int i=0; i<n_state; i++) rbfif[i] = new double[n_time_instance];
		for (int i=0; i<n_state; i++) for (int f=0; f<n_time_instance; f++) rbfif[i][f] = NEG_INFINITY;
	}

	for (int i=0; i<n_state; i++) for (int f=0; f<n_time_instance; f++) {
		double temp = NEG_INFINITY;
		for (int j=0; j<n_state; j++) 
			temp = max(temp, rbfs[i][j][0][f]);
		rbfif[i][f] = temp;
	}

	// set rbfjs
	if (rbfjs == NULL) { // allocate the memory for rbf_{j}[s,H)
		rbfjs = new double*[n_state];
		for (int j=0; j<n_state; j++) rbfjs[j] = new double[n_time_instance];
		for (int j=0; j<n_state; j++) for (int s=0; s<n_time_instance; s++) rbfjs[j][s] = NEG_INFINITY;
	}

	for (int j=0; j<n_state; j++) for (int s=0; s<n_time_instance; s++) {
		double temp = NEG_INFINITY;
		for (int i=0; i<n_state; i++) 
			temp = max(temp, rbfs[i][j][s][n_time_instance-1]);
		rbfjs[j][s] = temp;
	}

	// generate execution request matrix
	if (exec_req_matrix == NULL) { // allocate memory
		exec_req_matrix = new double*[n_state];
		for (int i=0; i<n_state; i++) exec_req_matrix[i] = new double[n_state];
	}

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		exec_req_matrix[i][j] = rbfs[i][j][0][n_time_instance-1];
	}
	exec_req_matrix_power[1] = exec_req_matrix;

	// delete rbfs
	if (rbfs != NULL) {
		for (int i=0; i<n_state; i++) {
			for (int j=0; j<n_state; j++) {
				for (int s=0; s< n_time_instance; s++) {
					delete[] rbfs[i][j][s];
				}
				delete[] rbfs[i][j];
			}
			delete[] rbfs[i];
		}
		delete[] rbfs;
	}
}

void Stateflow::generate_ibfs() {
	// Sort the transitions by the priorities
	// sort(trans.begin(),trans.end(),less_compare_trans);
	// allocate memory for rbf_{i,j}[s,f)
	rbfs = new double***[n_state];
	for (int i=0; i<n_state; i++) rbfs[i] = new double**[n_state];
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) rbfs[i][j] = new double*[n_time_instance];
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) for (int k=0; k<n_time_instance; k++) rbfs[i][j][k] = new double[n_time_instance];

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++)
		for (int s=0; s<n_time_instance; s++) for (int f=0; f<n_time_instance; f++) {
			if (i==j && s<=f) rbfs[i][j][s][f] = 0;
			else rbfs[i][j][s][f] = NEG_INFINITY;
		}
	
	// Dynamic Programming
	for (int len=1; len<n_time_instance; len++) {
		
		for (int s=0; s<n_time_instance-len; s++) {
			int f = s+len; // index of finish time
			int prev_ft = index_time[f-1]; // previous finish time

			//initialize
			for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
				if (len==1) {
					if (i==j) rbfs[i][j][s][f]=0;
				} else rbfs[i][j][s][f] = rbfs[i][j][s][f-1];
			}

			// calculate rbfs[i][j][s][f] based on the assumption we have known all the rbf[x][y][s][s..f-1].

			for (vector<Transition*>::iterator t_iter = trans.begin(); t_iter != trans.end(); t_iter++) {
				Transition* tran = *t_iter;
				int srcIndx = state_index[tran->src];
				int snkIndx = state_index[tran->snk];

				if(prev_ft%tran->period == 0) { // this transition might happen
					for (vector<State*>::iterator s_iter = states.begin(); s_iter != states.end(); s_iter++) {
						State* state = *s_iter;
						int sIndx = state_index[state];

						if (len == 1 && sIndx == srcIndx)
							rbfs[sIndx][snkIndx][s][f] = max(rbfs[sIndx][snkIndx][s][f], (double)tran->wcet);
						else if (len >1)
							rbfs[sIndx][snkIndx][s][f] = max(rbfs[sIndx][snkIndx][s][f], rbfs[sIndx][srcIndx][s][f-1]+tran->wcet);
					}
				}
			}
		}
	}

	// set rbfsf
	if (rbfsf == NULL) { // allocate the memory for rbf[s,f)
		rbfsf = new double*[n_time_instance];
		for (int s=0; s<n_time_instance; s++) rbfsf[s] = new double[n_time_instance]; 
		for (int s=0; s<n_time_instance; s++) for (int f=0; f<n_time_instance; f++) {
			if ( s <= f) rbfsf[s][f] = 0.0;
			else rbfsf[s][f] = NEG_INFINITY;
		}
	}

	for (int s=0; s<n_time_instance; s++) for (int f=0; f<n_time_instance; f++) {
		double temp = NEG_INFINITY;
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
			temp = max(temp, rbfs[i][j][s][f]);
		}
		rbfsf[s][f] = temp;
	}


	// set ibfsf
	for (int s=0; s<n_time_instance; s++) {
		map<int,double> map_ibfs;
		int start = index_time[s];
		
		for (int f=s; f<n_time_instance; f++) {
			if (s == f) continue;
			int finish = index_time[f];

			map<int,double> temp = map_calculate_ibf_within_one_hyperperiod(start,finish);
			for (map<int,double>::iterator iter = temp.begin(); iter != temp.end(); iter++) {
				map_ibfs[iter->first] = iter->second;
			}
		}

		map_ibfs[start] = 0;

		// output map_ibfs
		if (false) {
			for (map<int,double>::iterator iter = map_ibfs.begin(); iter != map_ibfs.end(); iter++) {
				cout << "ibf[" << start << "," << iter->first << ")=" << iter->second << endl;
			}
		}

		ibfsf[start] = map_ibfs;
	}

	// set rbfif
	if (rbfif == NULL) { // allocate the memory for rbf_{i}[0,f)
		rbfif = new double*[n_state];
		for (int i=0; i<n_state; i++) rbfif[i] = new double[n_time_instance];
		for (int i=0; i<n_state; i++) for (int f=0; f<n_time_instance; f++) rbfif[i][f] = NEG_INFINITY;
	}

	for (int i=0; i<n_state; i++) for (int f=0; f<n_time_instance; f++) {
		double temp = NEG_INFINITY;
		for (int j=0; j<n_state; j++) 
			temp = max(temp, rbfs[i][j][0][f]);
		rbfif[i][f] = temp;
	}

	// set ibfif
	for (int i=0; i<n_state; i++) {
		map<int,double> ret;
		int start = 0;

		for (int f=1; f<n_time_instance; f++) {
			int finish = index_time[f];

			for (int j=0; j<n_state; j++) {
				map<int,double> temp = map_calculate_ibf_within_one_hyperperiod(i,j,start,finish);

				if (ret.empty()) {
					ret=temp;
					continue;
				}

				for (map<int,double>::iterator iter = temp.begin(); iter != temp.end(); iter++) {

					int x = iter->first;
					int y = iter->second;

					bool dominated = false;
					bool dominating = false;
					set<int> remove;

					for (map<int,double>::iterator iter1 = ret.begin(); iter1 != ret.end(); iter1++) {
						int x1 = iter1->first;
						int y1 = iter1->second;

						if (y > y1 && y-y1 >= x-x1) {
							dominating = true;
							remove.insert(x1);
						}

						if (y1 >= y && y1-y >= x1-x) {
							dominated = true;
						}
					}

					if(dominating) {
						// delete the elements dominated
						for (set<int>::iterator iter1 = remove.begin(); iter1 != remove.end(); iter1++) {
							ret.erase(*iter1);
						}
					}

					if (!dominated)
						ret[x] = y;
				}
			}
		}

		map<int,double> map_ibfi;
		for (map<int,double>::iterator iter = ret.begin(); iter != ret.end(); iter++) {
			double y = iter->second;
			map<int,double>::iterator next_Iter = iter;
			next_Iter++;
			int x;
			if (next_Iter == ret.end()) {
				x = hyperperiod;
			}
			else {
				x = next_Iter->first - ceil(next_Iter->second - iter->second);
			}

			map_ibfi[x] = y;
		}

		map_ibfi[start] = 0;

		// output map_ibfs
		if (false) {
			for (map<int,double>::iterator iter = map_ibfi.begin(); iter != map_ibfi.end(); iter++) {
				cout << "ibf_{"<<i<<"}[" << start << "," << iter->first << ")=" << iter->second << endl;
			}
		}

		ibfif[i] = map_ibfi;
	}
	
	// set ibfjs, it should equal rbfjs
	if (ibfjs == NULL) { // allocate the memory for ibf_{j}[s,H)
		ibfjs = new double*[n_state];
		for (int j=0; j<n_state; j++) ibfjs[j] = new double[n_time_instance];
		for (int j=0; j<n_state; j++) for (int s=0; s<n_time_instance; s++) ibfjs[j][s] = NEG_INFINITY;
	}

	for (int j=0; j<n_state; j++) for (int s=0; s<n_time_instance; s++) {
		double temp = NEG_INFINITY;
		for (int i=0; i<n_state; i++) 
			temp = max(temp, rbfs[i][j][s][n_time_instance-1]);
		ibfjs[j][s] = temp;
	}

	// generate execution request matrix
	if (exec_req_matrix == NULL) { // allocate memory
		exec_req_matrix = new double*[n_state];
		for (int i=0; i<n_state; i++) exec_req_matrix[i] = new double[n_state];
	}

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		exec_req_matrix[i][j] = rbfs[i][j][0][n_time_instance-1];
	}
	exec_req_matrix_power[1] = exec_req_matrix;

	// delete rbfs
	if (rbfs != NULL) {
		for (int i=0; i<n_state; i++) {
			for (int j=0; j<n_state; j++) {
				for (int s=0; s< n_time_instance; s++) {
					delete[] rbfs[i][j][s];
				}
				delete[] rbfs[i][j];
			}
			delete[] rbfs[i];
		}
		delete[] rbfs;
	}
}

void Stateflow::generate_exec_req_matrix() {
	// We have implemented this function in generate_rbfs or generate_ibfs
}

void Stateflow::calculate_exec_req_matrix_power(int tf) {
	Timer timer;
	/// Calculate all the execution request matrix from 2 to tf
	for (int t=2; t<tf; t++) {
		if (exec_req_matrix_power.find(t) != exec_req_matrix_power.end())
			continue;
		double** temp;
		if (t>gdef+gper && isIrred && exec_req_matrix_power.find(t-gper) != exec_req_matrix_power.end()) {
			//cout<<"a-ha!"<<endl;
			//timer.start();
			temp = MaxPlusAlgebra::periodicly_calculate_maxplus_matrix_power(exec_req_matrix_power[t-gper],n_state,lfac*hyperperiod,gper);
			//timer.end();
			//double t0 = timer.getTime();

			/*
			timer.start();
			//double** temp1 = MaxPlusAlgebra::multiply_maxplus_matrix(exec_req_matrix_power[t-1],exec_req_matrix,n_state);
			double** temp1 = MaxPlusAlgebra::power_maxplus_matrix(exec_req_matrix,n_state, t);
			timer.end();
			double t1 = timer.getTime();

			for (int i=0; i<n_state; i++) delete[] temp1[i];
			delete[] temp1;

			//cout<<t0<<"\t"<<t1<<"\t"<<t1-t0<<endl;

			tDiff += t1-t0;
			tCalWithoutPeriodicity += t1;
			tCalWithPeriodicity += t0;
			*/
		}
		else
			temp = MaxPlusAlgebra::multiply_maxplus_matrix(exec_req_matrix_power[t-1],exec_req_matrix,n_state);
		exec_req_matrix_power[t] = temp;
	}
}

/*
void Stateflow::calculate_exec_req_matrix_power_DP_with_periodicity(int tf) {

}

void Stateflow::calculate_exec_req_matrix_power_DP_without_periodicity(int tf) {
	int nSize = states.size();
	int nPeriod = tf;

	

	double** preNodeReqBoundFunc = new double*[nSize*nSize];
	for (int i=0; i<nSize*nSize; i++) preNodeReqBoundFunc[i] = new double[nPeriod+1];

	for (int i=0; i<nSize*nSize; i++) {
		int x = i/nSize;
		int y = i%nSize;
		preNodeReqBoundFunc[i][0] = exec_req_matrix[x][y];
	}

	for (int i=1; i<=nPeriod; i++) {
		int release = i;

		for (vector<Node*>::iterator iter = node_vec.begin(); iter != node_vec.end(); iter++) {
			Node* node = *iter;
			int nodeIndex = node_to_index[node];
			
			double nodeRbf = node->wcet;
			for (list<Edge*>::iterator iter2 = node->in.begin(); iter2 != node->in.end(); iter2++) {
				Edge* edge = *iter2;
				
				int srcRelease = (release-edge->separationTime)/this->pGCD;
				if (srcRelease >= 0) {
					int srcIndex = node_to_index[edge->src_node];
					nodeRbf = max(nodeRbf, node->wcet + preNodeReqBoundFunc[srcIndex][srcRelease]);
				}
			}
			preNodeReqBoundFunc[nodeIndex][i] = nodeRbf;
		}

		
	}

	for (int i=0; i<=nPeriod; i++) {
		int release = i*this->pGCD;
		if (i==0) {
			rbfs[release] = 0;
			continue;
		}

		double noderbf = 0;
		for (int j=0; j<nSize; j++) {
			noderbf = max(noderbf, preNodeReqBoundFunc[j][i-1]);
		}

		rbfs[release] = noderbf;
	}

	for (int i=0; i<nSize; i++) {
		//delete[] nodeReqBoundFunc[i];
		delete[] preNodeReqBoundFunc[i];
	}

	//delete[] nodeReqBoundFunc;
	delete[] preNodeReqBoundFunc;
}
*/
void Stateflow::calculate_generialized_factor() {
	gfac = Utility::creat_matrix(n_state,n_state);
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) gfac[i][j] = lfac;
}

void Stateflow::calculate_generialized_period() {
	if (!isIrred) return;
	gper = MaxPlusAlgebra::calculate_linear_period(exec_req_matrix,n_state,lfac*hyperperiod,false);
}

void Stateflow::calculate_generialized_defect() {
	if (!isIrred) 
		gdef = 100*scale;
	else
		gdef = MaxPlusAlgebra::calculate_linear_defect(exec_req_matrix_power, n_state, lfac*hyperperiod, gper, tf0);
}

void Stateflow::check_irreducible() {
	//this->lfac = precise_digraph->linear_factor;
	this->isIrred = simple_digraph->strongly_connected;
}

void Stateflow::calculate_linear_factor() {
	if (isIrred)
		lfac = GraphAlgorithms::calculate_maximum_cycle_mean(exec_req_matrix, n_state)/hyperperiod;
	else {
		exec_digraph = new GeneralDirectedGraph(n_state, exec_req_matrix);
		//write_graphviz(cout);
		//exec_digraph->write_graphviz(cout);
		exec_digraph->generate_strongly_connected_components();
		if (exec_digraph->strongly_connected) isIrred = true;
		exec_digraph->calculate_untilization();
		lfac = exec_digraph->util/hyperperiod;
		delete exec_digraph;
	}
}

void Stateflow::scale_wcet(double factor) {
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* t = *iter;
		int wcet = (int)(factor*t->wcet);
		if (wcet == 0) wcet = rand()%5+1; // force the wcet not to be 0
		t->wcet = wcet;
	}
}

void Stateflow::calculate_csum() {
	int _csum = 0;
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		//int temp = tran->wcet*(hyperperiod/tran->period);
		_csum += tran->wcet*(hyperperiod/tran->period);
	}
	csum = _csum;
}

void Stateflow::calculate_tf0(Stateflow** stateflows, int i) {
	double util = 0;
	double sum = 0;
	for (int k=0; k<=i; k++) {
		Stateflow* sf = stateflows[k];
		util += sf->lfac;

		sum += sf->csum*2;
		if (k==i) sum -= sf->csum;
	}

	if (util >= 1) 
		tf0 = POS_INFINITY;
	else 
		tf0 = sum/(1-util);
}

/// Generate a simple digraph for stateflow
void Stateflow::generate_simple_digraph() {
	if (simple_digraph != NULL) return;

	simple_digraph = GraphAlgorithms::generate_simple_digraph(this);
	// prepare for calculating linear upper bounds
	// simple_digraph->prepare_digraph();
}

/// Generate a precise digraph for stateflow
void Stateflow::generate_precise_digraph() {
	if (precise_digraph != NULL) return;
	precise_digraph = GraphAlgorithms::generate_precise_digraph(this);

	// prepare for calculating linear upper bounds
	// precise_digraph->prepare_digraph();

	//lfac = precise_digraph->linear_factor;
}

/// Generate a strongly-connected precise digraph for stateflow
void Stateflow::generate_strongly_connected_precise_digraph() {
	if (sc_precise_digraph != NULL) return;
	sc_precise_digraph = GraphAlgorithms::generate_strongly_connected_precise_digraph(this);
}

/// Generate a reachability digraph for stateflow
void Stateflow::generate_reachability_digraph() {
	if (reachability_digraph != NULL) return;
	reachability_digraph = GraphAlgorithms::generate_reachability_digraph(this);
}

/// Calculate the tigher linear upper bounds
/// true: operating on precise digraph
/// false: operating on simple digraph
void Stateflow::calculate_linear_upper_bounds(bool precise) {
	if(precise) {
		// generate_precise_digraph();
		precise_digraph->calculate_all_gcd();
		precise_digraph->calculate_period_gcd();
		precise_digraph->generate_strongly_connected_components();
		precise_digraph->check_strongly_connected();
		precise_digraph->calculate_linear_factor();
		precise_digraph->calculate_linear_upper_bounds();

		lfac = precise_digraph->linear_factor;
		crbf = precise_digraph->c_rbf;
		cibf = precise_digraph->c_ibf;
		cdbf = precise_digraph->c_dbf;
	} else {
		//  generate_simple_digraph();
		simple_digraph->calculate_all_gcd();
		simple_digraph->calculate_period_gcd();
		simple_digraph->generate_strongly_connected_components();
		simple_digraph->check_strongly_connected();
		simple_digraph->calculate_linear_factor();
		simple_digraph->calculate_linear_upper_bounds();

		lfac = simple_digraph->linear_factor;
		crbf = simple_digraph->c_rbf;
		cdbf = simple_digraph->c_dbf;
		cibf = simple_digraph->c_ibf;
	}
}

void Stateflow::calculate_tf1(Stateflow** stateflows, int i) {
	double util = stateflows[i]->lfac;
	double sum = stateflows[i]->cdbf;
	for (int k=0; k<i; k++) {
		Stateflow* sf = stateflows[k];
		util += sf->lfac;

		sum += sf->crbf;
	}

	if (util >= 1) 
		tf1 = POS_INFINITY;
	else
		tf1 = sum/(1-util);
}

void Stateflow::calculate_tf2(Stateflow** stateflows, int i) {
	double util = stateflows[i]->lfac;
	double sum = stateflows[i]->cdbf;
	for (int k=0; k<i; k++) {
		Stateflow* sf = stateflows[k];
		util += sf->lfac;

		sum += sf->cibf;
	}

	if (util >= 1) 
		tf2 = POS_INFINITY;
	else
		tf2 = sum/(1-util);
}

/// Static offset
double Stateflow::get_rbf(int start, int finish) { // return rbf[s,f)
	if (start > finish) return NEG_INFINITY;
	if (start == finish) return 0;

	int s = start%hyperperiod;
	int sprim = *rbf_time_instances.lower_bound(s);
	int ns = start/hyperperiod;
	s = sprim;

	int f = finish%hyperperiod;
	int fprim = *rbf_time_instances.lower_bound(f);
	int nf = finish/hyperperiod;
	int n = nf - ns;
	f = n*hyperperiod+fprim;

	for (vector<StartFinishTime*>::iterator iter = rbf_vec.begin(); iter != rbf_vec.end(); iter++) {
		StartFinishTime* sft = *iter;
		if (sft->start == s && sft->finish == f)
			return sft->value;
	}
	
	double rt = 0;

	if (f<=hyperperiod) // s and f are in the same hyperperiod
		rt = calculate_rbf_within_one_hyperperiod(s,f);
	else 
		rt = calculate_rbf_within_multiple_hyperperiods(s,f);

	StartFinishTime* sft = new StartFinishTime(s,f,rt);
	rbf_vec.push_back(sft);

	return rt;
}

double Stateflow::calculate_rbf_within_one_hyperperiod(int start, int finish) {
	/*
	double rt = NEG_INFINITY;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		rt = max(rt, calculate_rbf_within_one_hyperperiod(i,j,start,finish));
	}
	return rt;
	*/

	int sIndx = time_index[start];
	int fIndx = time_index[finish];

	return rbfsf[sIndx][fIndx];
}

double Stateflow::calculate_rbf_within_one_hyperperiod(int i, int j, int start, int finish) {
	int sIndx = time_index[start];
	int fIndx = time_index[finish];

	return rbfs[i][j][sIndx][fIndx];
}

double Stateflow::calculate_rbf_within_multiple_hyperperiods(int start, int finish) {

	int n = finish/hyperperiod;
	int f = finish - n*hyperperiod;

	//double** prev = Utility::creat_matrix(n_state,n_state);
	double** midd = Utility::creat_matrix(n_state,n_state);
	//double** post = Utility::creat_matrix(n_state,n_state);

	double* prev = new double[n_state];
	double* post = new double[n_state];

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		//prev[i][j] = 0;
		midd[i][j] = 0;
		//post[i][j] = 0;
	}
	
	int sIndx = time_index[start];
	int fIndx = time_index[f];

	/*
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		prev[i][j] = calculate_rbf_within_one_hyperperiod(i,j,start,hyperperiod);
	}
	*/
	for (int j=0; j<n_state; j++)
		prev[j] = rbfjs[j][sIndx];

	for (int i=0; i<n_state; i++)
		post[i] = rbfif[i][fIndx];

	if (n>1) {
		if (exec_req_matrix_power.find(n)==exec_req_matrix_power.end()) 
			calculate_exec_req_matrix_power(n);
		
		for (int i=0; i<n_state; i++) for (int j=0; j<n_state;j++)
			midd[i][j] = exec_req_matrix_power[n-1][i][j];
	}

	double rbf = 0;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		double rbfij = 0;
		if (n!=1) {
			for (int k=0; k<n_state; k++) for (int l=0; l<n_state; l++) {
				rbfij = max(rbfij, prev[k]+midd[k][l]+post[l]);
			}
		}
		else {// n==1, we should not use the middle matrix in case of the errors of states
			for (int k=0; k<n_state; k++) {
				rbfij = max(rbfij, prev[k]+post[k]);
			}
		}
		rbf = max(rbf, rbfij);
	}

	// release prev, midd and post matrices
	for (int i=0; i<n_state; i++) {
		//delete[] prev[i];
		delete[] midd[i];
		//delete[] post[i];
	}
	delete[] prev;
	delete[] midd;
	delete[] post;

	return rbf;
}

double Stateflow::get_ibf(int start, int finish) { // return ibf[s,f)
	if (start > finish) return NEG_INFINITY;
	if (start == finish) return 0;

	int s = start%hyperperiod;
	int sprim = *rbf_time_instances.lower_bound(s);
	int ns = start/hyperperiod;
	s = sprim;

	int f = finish%hyperperiod;
	int nf = finish/hyperperiod;
	int n = nf - ns;
	f = n*hyperperiod+f;

	for (vector<StartFinishTime*>::iterator iter = ibf_vec.begin(); iter != ibf_vec.end(); iter++) {
		StartFinishTime* sft = *iter;
		if (sft->start == s && sft->finish == f)
			return sft->value;
	}
	
	double rt = 0;

	if (f<=hyperperiod) { // s and f are in the same hyperperiod
		map<int,double> map_ibfs = ibfsf[s];
		map<int,double>::const_iterator lb_iter = map_ibfs.lower_bound(f);
		if (lb_iter == map_ibfs.end()) {
			cerr << "There is no f=" << f << "\t in map_ibfs" << endl;
			// output map_ibfs
			for (map<int,double>::iterator iter = map_ibfs.begin(); iter != map_ibfs.end(); iter++) {
				cerr << "ibf[" << start << "," << iter->first << ")=" << iter->second << endl;
			}
			exit(EXIT_FAILURE);
		}

		rt = lb_iter->second;
	}
	else 
		rt = calculate_ibf_within_multiple_hyperperiods(s,f);

	StartFinishTime* sft = new StartFinishTime(s,f,rt);
	ibf_vec.push_back(sft);

	return rt;
}

double Stateflow::calculate_ibf_within_one_hyperperiod(int start, int finish) {
	double rt = NEG_INFINITY;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		rt = max(rt, calculate_ibf_within_one_hyperperiod(i,j,start,finish));
	}
	return rt;
}

double Stateflow::calculate_ibf_within_one_hyperperiod(int i, int j, int start, int finish) {
	int sIndx = time_index[start];
	int fIndx = 0;

	int fh = *rbf_time_instances.lower_bound(finish);
	if (fh == finish) {
		fIndx = time_index[finish];
		return rbfs[i][j][sIndx][fIndx];
	}
	
	int fl = *(--rbf_time_instances.lower_bound(finish));
	int fIndx0 = time_index[fl];

	double ibfij = rbfs[i][j][sIndx][fIndx0];
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* tran = *iter;
		int src = state_index[tran->src];
		int snk = state_index[tran->snk];
		if (snk != j) continue;
		if (fl%tran->period==0) {
			double dmin = min(tran->wcet, finish-fl);
			double rbfij = rbfs[i][src][sIndx][fIndx0];
			if (rbfij == NEG_INFINITY) continue;
			ibfij = max(ibfij, rbfij+dmin);
		}
	}

	return ibfij;
}

map<int,double> Stateflow::map_calculate_ibf_within_one_hyperperiod(int start, int finish) {
	map<int,double> m_ibfij;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		map<int,double> temp = map_calculate_ibf_within_one_hyperperiod(i,j,start,finish);
		if (i == 0 && j == 0) {
			m_ibfij=temp;
			continue;
		}

		for (map<int,double>::iterator iter = temp.begin(); iter != temp.end(); iter++) {

			int x = iter->first;
			int y = iter->second;

			bool dominated = false;
			bool dominating = false;
			set<int> remove;

			for (map<int,double>::iterator iter1 = m_ibfij.begin(); iter1 != m_ibfij.end(); iter1++) {
				int x1 = iter1->first;
				int y1 = iter1->second;

				if (y > y1 && y-y1 >= x-x1) {
					dominating = true;
					remove.insert(x1);
				}

				if (y1 >= y && y1-y >= x1-x) {
					dominated = true;
				}
			}

			if(dominating) {
				// delete the elements dominated
				for (set<int>::iterator iter1 = remove.begin(); iter1 != remove.end(); iter1++) {
					m_ibfij.erase(*iter1);
				}
			}

			if (!dominated)
				m_ibfij[x] = y;
		}
	}
	
	map<int,double> ret;
	for (map<int,double>::iterator iter = m_ibfij.begin(); iter != m_ibfij.end(); iter++) {
		double y = iter->second;
		map<int,double>::iterator next_Iter = iter;
		next_Iter++;
		int x;
		if (next_Iter == m_ibfij.end()) {
			x = finish;
		}
		else x = next_Iter->first;

		ret[x] = y;
	}

	/*
	for (map<int,double>::iterator iter = ret.begin(); iter != ret.end(); iter++) {
		cout << "map_calculate_ibf_within_one_hyperperiod=>" << "ibf[" << start << "," << iter->first << "]=" << iter->second << endl;
	}
	*/

	return ret;
}

map<int,double> Stateflow::map_calculate_ibf_within_one_hyperperiod(int i, int j, int start, int finish) {
	// start and finish should belong the set rbf_time_instances
	int sIndx = time_index[start];
	int fIndx = time_index[finish];

	State* sj = index_state[j];

	//cout << "finish=" << finish << endl;
	
	int fl = *(--rbf_time_instances.lower_bound(finish));
	//cout << "f1" << fl << endl;
	int fIndx0 = time_index[fl];

	map<int,double> ibfRec;

	for (int k=0; k<n_state; k++) {
		double ibfik = rbfs[i][k][sIndx][fIndx0];
		if (ibfik == NEG_INFINITY) continue;
		State* sk = index_state[k];
		int wcet = 0;// sk->getMaximalWCET(sj,fl);

		for (list<Transition*>::iterator iter = sk->out.begin(); iter != sk->out.end(); iter++) {
			if ((*iter)->snk != sj) continue;
			if (fl%(*iter)->period != 0) continue;
			wcet = max(wcet, (*iter)->wcet);
		}

		if (ibfRec.find(fl) == ibfRec.end()) {
			ibfRec[fl] = ibfik;
		}
		else {
			ibfRec[fl] = max(ibfRec[fl],ibfik);
		}

		int x = fl + wcet;
		int y = ibfik+wcet;

		bool dominated = false;
		bool dominating = false;
		set<int> remove;

		for (map<int,double>::iterator iter = ibfRec.begin(); iter != ibfRec.end(); iter++) {
			int x1 = iter->first;
			int y1 = iter->second;

			if (y > y1 && y-y1 >= x-x1) {
				dominating = true;
				remove.insert(x1);
			}

			if (y1 >= y && y1-y >= x1-x) {
				dominated = true;
			}
		}

		if(dominating) {
			// delete the elements dominated
			for (set<int>::iterator iter = remove.begin(); iter != remove.end(); iter++) {
				ibfRec.erase(*iter);
			}
		}

		if (!dominated)
			ibfRec[x] = y;
	}	

	//ibfRec.erase(start);
	/*
	for (map<int,double>::iterator iter = ibfRec.begin(); iter != ibfRec.end(); iter++) {
		cout << "map_calculate_ibf_within_one_hyperperiod=>" << "ibf_{" << i << "," << j << "}[" << start << "," << iter->first << "]=" << iter->second << endl;
	}
	*/

	return ibfRec;
}

double Stateflow::calculate_ibf_within_multiple_hyperperiods(int start, int finish) {

	int n = finish/hyperperiod;
	int f = finish - n*hyperperiod;

	//double** prev = Utility::creat_matrix(n_state,n_state);
	double** midd = Utility::creat_matrix(n_state,n_state);
	//double** post = Utility::creat_matrix(n_state,n_state);

	double* prev = new double[n_state];
	double* post = new double[n_state];

	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		//prev[i][j] = 0;
		midd[i][j] = 0;
		//post[i][j] = 0;
	}

	int sIndx = time_index[start];
	
	for (int j=0; j<n_state; j++)
		prev[j] = ibfjs[j][sIndx];

	for (int i=0; i<n_state; i++) {
		map<int,double> ibfi = ibfif[i];
		map<int,double>::const_iterator ibfi_iter = ibfi.lower_bound(f);
		if (ibfi_iter == ibfi.end()) {
			cerr << "No f=" << f << " in ibf_{" << i << "}" << endl;
			exit(EXIT_FAILURE);
		}
		post[i] = ibfi_iter->second;
	}
	

	if (n>1) {
		if (exec_req_matrix_power.find(n)==exec_req_matrix_power.end()) 
			calculate_exec_req_matrix_power(n);

		for (int i=0; i<n_state; i++) for (int j=0; j<n_state;j++)
			midd[i][j] = exec_req_matrix_power[n-1][i][j];
	}

	double ibf = 0;
	for (int i=0; i<n_state; i++) for (int j=0; j<n_state; j++) {
		double ibfij = 0;
		if (n!=1) {
			for (int k=0; k<n_state; k++) for (int l=0; l<n_state; l++) {
				ibfij = max(ibfij, prev[k]+midd[k][l]+post[l]);
			}
		}
		else {// n==1, we should not use the middle matrix in case of the errors of states
			for (int k=0; k<n_state; k++) {
				ibfij = max(ibfij, prev[k]+post[k]);
			}
		}
		ibf = max(ibf, ibfij);
	}

	// release prev, midd and post matrices
	for (int i=0; i<n_state; i++) {
		//delete[] prev[i];
		delete[] midd[i];
		//delete[] post[i];
	}
	delete[] prev;
	delete[] midd;
	delete[] post;

	return ibf;
}

/// Note: dbf[s,f] = rbf[s,lower(f))
double Stateflow::get_dbf(int start, int finish) { // return dbf[s,f)
	if (start > finish) return NEG_INFINITY;
	if (start == finish) return 0;

	int s = start%hyperperiod;
	int sprim = *rbf_time_instances.lower_bound(s);
	int ns = start/hyperperiod;
	s = sprim;

	int f = finish%hyperperiod;
	int hf = *rbf_time_instances.lower_bound(f);
	int hl = 0;
	if (rbf_time_instances.find(f) != rbf_time_instances.end())
		hl = f;
	else {
		if (f!=0) hl = *(--rbf_time_instances.lower_bound(f));
	}
	int nf = finish/hyperperiod;

	int n = nf - ns;

	/*
	if (f == hf) {
		f = n*hyperperiod+hf;
	} 
	else if (f-hl >= gcd) {
		f = n*hyperperiod+hf;
	} 
	else {
		f = n*hyperperiod+hl;
	}
	*/

	f = n*hyperperiod+hl;
	//cout << finish << "\t" << f << "\t" << hf << "\t" << hl << endl;
	
	if (s > f) return NEG_INFINITY;
	if (s == f) return 0;

	return get_rbf(s,f);

}

/// Arbitrary offset
double Stateflow::get_rbf(int t) { // return rbf(t)
	if (t==0) return 0;

	if (map_rbf.find(t) != map_rbf.end())
		return map_rbf[t];

	double rbf = 0;
	for (int i=0; i<n_time_instance; i++) {
		int s = index_time[i];
		int f = s+t;
		double rbfsf = get_rbf(s,f);
		rbf = max(rbf, rbfsf);
	}

	map_rbf[t] = rbf;

	return rbf;
}

double Stateflow::get_ibf(int t) { // return ibf(t)
	if (t==0) return 0;
	
	if (map_ibf.find(t) != map_ibf.end())
		return map_ibf[t];

	double ibf = 0;
	for (int i=0; i<n_time_instance; i++) {
		int s = index_time[i];
		int f = s+t;
		double ibfsf = get_ibf(s,f);
		ibf = max(ibf, ibfsf);
	}

	map_ibf[t] = ibf;

	return ibf;
}

double Stateflow::get_dbf(int t) { // return dbf(t)
	if (t==0) return 0;

	if (map_dbf.find(t) != map_dbf.end())
		return map_dbf[t];

	double dbf = 0;
	for (int i=0; i<n_time_instance; i++) {
		int s = index_time[i];
		int f = s+t;
		double dbfsf = get_dbf(s,f);
		dbf = max(dbf, dbfsf);
	}

	map_dbf[t] = dbf;

	return dbf;
}

void Stateflow::write_graphviz(ostream& out) {
	out << "digraph G {" <<endl;
	for (vector<State*>::iterator iter = states.begin(); iter != states.end(); iter++) {
		State* state = *iter;
		out << state->name << " [label=\" " << state->name << " \"]" << endl;
	}
	for (vector<Transition*>::iterator iter = trans.begin(); iter != trans.end(); iter++) {
		Transition* transition = *iter;
		out << transition->src->name << " -> " << transition->snk->name 
			<< " [label=\" " << transition->wcet << " / " << transition->period << " \"]" << endl;
	}
	out << "}" <<endl;
}


