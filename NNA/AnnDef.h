#ifndef ANN_HEADER
#define ANN_HEADER
#include <iostream>
#include <Eigen/dense>
#include <random> 
#include <algorithm>
#include "AnnImp.cpp"
#define NPOP 100
#define MAX  140
//float L1 = 0.05, U1 = 0.20, L2 = 0.25, U2 = 1.3, L3 = 10.00, U3 = 12.00;//tension
//float L1 = 0.0, U1 = 99.0, L2 = 0.0, U2 = 99.0, L3 = 10.00, U3 = 100.0, L4 = 10.00, U4 = 100.0; //pressure vessel
//float L1 = 47.0, U1 = 51.0, L2 = 75.0, U2 = 85.0, L3 = 0.99, U3 = 5.0, L4 = 0.99, U4 = 5.0;//i-beam 
float L1 = 0.1, U1 = 2.0, L2 = 0.1, U2 = 10.0, L3, U3, L4, U4;//i-beam 

#define nvars 4

namespace NNA {
	static double f(Eigen::Matrix<double, 1, nvars>);
	static void initialization(Eigen::Matrix<double, 1, NPOP>& ww, Eigen::Matrix<double, NPOP, NPOP>& w,
		Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& cost);
	static void CreateInitialPopulation(Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& cost,
		Eigen::Matrix <double, 1, 2>& value_index);
	static void generateWeight(Eigen::Matrix<double, NPOP, NPOP>& w, Eigen::Matrix<double, 1, NPOP - 1>& t);
	static void setTarget(Eigen::Matrix<double, 1, nvars>& x_target, Eigen::Matrix<double, NPOP, nvars>& x_pattern,
		Eigen::Matrix <double, 1, 2>& value_index, Eigen::Matrix<double, NPOP, 1>& w_target, Eigen::Matrix<double, NPOP, NPOP>& w);
	static void Run(Eigen::Matrix<double, NPOP, nvars>& x_new, Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& w_target,
		Eigen::Matrix<double, NPOP, NPOP>& w, Eigen::Matrix<double, NPOP, 1>& cost, Eigen::Matrix<double, 1, nvars>& x_target);
};
#endif