#include "Utility.h"

double Utility::EPSILON = 0.000001;

extern double POS_INFINITY;
extern double NEG_INFINITY;

// return the greatest common divisor for two integers
int Utility::math_gcd(int a, int b) {
		if (b==0) return a;
		return math_gcd(b, a%b);
}

int Utility::math_lcm(int a, int b) {
		if (b==0) return a;
		return (a*(b/math_gcd(a,b)));
}

std::string Utility::int_to_string(int a) {
		char temp[10];
		sprintf(temp,"%d",a);
		return temp;
}

double** Utility::creat_matrix(int nrow, int ncol) {
	double** A = new double*[nrow];
	for (int i=0; i<nrow; i++) A[i] = new double[ncol];
	return A;
}

void Utility::output_matrix(double** A, int nrow, int ncol) {
	for (int i=0; i<nrow; i++) {
		for (int j=0; j<ncol; j++) {
			std::cout<< A[i][j];
			if (j==ncol-1) std::cout<<std::endl;
			else std::cout<<"\t";
		}
	}
}

void Utility::output_latex_matrix(double** A, int nrow, int ncol) {
	for (int i=0; i<nrow; i++) {
		std::cout << "c";
	}
	std::cout << std::endl;

	for (int i=0; i<nrow; i++) {
		for (int j=0; j<ncol; j++) {
			if (A[i][j] == NEG_INFINITY) std::cout << "-\\infty";
			else std::cout<< 0.01*A[i][j];
			if (j==ncol-1) std::cout << "\\\\" <<std::endl;
			else std::cout<<"&";
		}
	}
}

double* Utility::uniformly_distributed(int n, double tUtil) {
	double* ret = new double[n];

	double sum = tUtil;
	for (int i=0; i<n-1; i++) {
		double nextsum = sum*pow((double)rand()/RAND_MAX, 1.0/(n-i-1));
		ret[i] = sum-nextsum;
		sum = nextsum;
	}

	ret[n-1] = sum;
	return ret;
}

double* Utility::uniformly_distributed2(int n, double tUtil) {
	double* ret = new double[n];

	double sum = tUtil;
	for (int i=0; i<n-1; i++) {
		double nextsum = sum*pow((double)rand()/RAND_MAX, 1.0/(n-i-1));
		if (nextsum < 0.001) nextsum = 0.001;
		nextsum = (int)(nextsum*1000);
		nextsum = nextsum/1000;
		ret[i] = sum-nextsum;
		sum = nextsum;
	}

	ret[n-1] = sum;
	return ret;
}

double Utility::max_uniformly_distributed(int n, double tUtil) {
	double* ret = new double[n];

	double sum = tUtil;
	for (int i=0; i<n-1; i++) {
		double nextsum = sum*pow((double)rand()/RAND_MAX, 1.0/(n-i-1));
		ret[i] = sum-nextsum;
		sum = nextsum;
	}

	ret[n-1] = sum;

	double max = 0;
	for (int i=0; i<n; i++) {
		if (max < ret[i]) max = ret[i];
	}

	delete[] ret;
	return max;
}

bool Utility::compare_two_matrices(double** A, double** B, int n) {
	for (int i=0; i<n; i++) for (int j=0; j<n; j++)
		if (abs(A[i][j]-B[i][j])>EPSILON) {
			if (false) {
				std::cout<<i<<","<<j<<std::endl;
				std::cout<<A[i][j]<<"..."<<B[i][j]<<std::endl;
			}
			return false;
		}
	return true;
}

void Utility::output_one_vector(ostream& out, string sVec, vector<double> vec) {
	out<<sVec<<"=[";
	typedef vector<double>::iterator Iter;
	for (Iter it = vec.begin(); it != vec.end(); it++) {
		out<<*it;
		if (it != --vec.end()) out<<",";
	}
	out<<"];"<<endl;
}

void Utility::output_one_vector(ostream& out, string sVec, vector<int> vec) {
	out<<sVec<<"=[";
	typedef vector<int>::iterator Iter;
	for (Iter it = vec.begin(); it != vec.end(); it++) {
		out<<*it;
		if (it != --vec.end()) out<<",";
	}
	out<<"];"<<endl;
}