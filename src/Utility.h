/* \file Utility.h
 *  This file implements some math-related functions.
 *  \author Chao Peng
 *  
 *  Changes
 *  ------
 *  22-Aug-2015 : Initial revision (CP)
 *
 */

#ifndef UTILITY_H
#define UTILITY_H

//#include <vld.h>

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <random>
#include <set>
#include <time.h>
#include <algorithm>
#include <functional> 

using namespace std;

class Utility {
public:
	static double EPSILON;

	// return the greatest common divisor for two integers
	static int math_gcd(int a, int b);
	static int math_lcm(int a, int b);
	static std::string int_to_string(int a);
	static double** creat_matrix(int nrow, int ncol);
	static void output_matrix(double** A, int nrow, int ncol);
	static void output_latex_matrix(double** A, int nrow, int ncol);
	static double* uniformly_distributed(int n, double tUtil);
	static double* uniformly_distributed2(int n, double tUtil);
	static double max_uniformly_distributed(int n, double tUtil);
	static bool compare_two_matrices(double** A, double** B, int n);
	static void output_one_vector(ostream& out, string sVec, vector<double> vec);
	static void output_one_vector(ostream& out, string sVec, vector<int> vec);
};

#endif