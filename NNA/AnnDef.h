#ifndef ANN_HEADER
#define ANN_HEADER
#include <iostream>
#include <Eigen/dense>
#include <random> 
#include <algorithm>
#include "AnnImp.cpp"
#define NPOP 100
#define MAX  140
#define LB 0.05 
#define UP 15.0
#define nvars 3

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