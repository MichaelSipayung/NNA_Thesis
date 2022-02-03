#ifndef ANN_HEADER
#define ANN_HEADER
#include <iostream>
#include <Eigen/dense>
#include <random> 
#include <algorithm>
#include "AnnImp.cpp"
#define NPOP 150
#define MAX  140
double tension[6] = {0.05,2.0,0.25,1.3,2.0,15.0};
double pressureVes[8] = { 1.0,99.0,1.0,99.0,10.00,200.0,10.0,200.0 };
double ibeam[8] = { 47.0,51.0,75.0,85.0,0.99,5.0,0.99,5.0 };
double weldeadBeam[4] = { 0.1,2.0,0.1,10.0 };
double LoSpeedRed[7] = {2.6,0.7,17.0,7.3,7.8,2.9,5.0};//speedReducer 
double UpSpeedRed[7] = {3.6,0.8,28.0,8.3,8.3,3.9,5.5 }; 
double cantBeam[2] = { 0.01,100.0 };
double currugatedBulkhead[4] = { 0.0,50.0,0.0,5.0 };
double tabColumn[2] = { 0.2 ,6.0 };
double gear[2] = { 12.0,60.0 };
double LoMulD[5] = { 60.0,90.0,1.5,600.0,2.0 };
double UpMulD[5] = { 80.0,110.0,3.0,1000.0,9.0 };
double bT[4] = { 1.5,2.0,2.5,3.0 }; //constrains for multiple disk clutch brake design 
double bF[41] = { 600.0 }; 
#define nvars 5
namespace NNA {
	static double f(Eigen::Matrix<double, 1, nvars>);
	static void initialization(Eigen::Matrix<double, 1, NPOP>& ww, Eigen::MatrixXd &w,
		Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& cost);
	static void CreateInitialPopulation(Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& cost,
		Eigen::Matrix <double, 1, 2>& value_index);
	static void generateWeight(Eigen::MatrixXd& w, Eigen::Matrix<double, 1, NPOP - 1>& t);
	static void setTarget(Eigen::Matrix<double, 1, nvars>& x_target, Eigen::Matrix<double, NPOP, nvars>& x_pattern,
		Eigen::Matrix <double, 1, 2>& value_index, Eigen::Matrix<double, NPOP, 1>& w_target, Eigen::MatrixXd& w);
	static void Run(Eigen::Matrix<double, NPOP, nvars>& x_new, Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& w_target,
		Eigen::MatrixXd &w, Eigen::Matrix<double, NPOP, 1>& cost, Eigen::Matrix<double, 1, nvars>& x_target);
};
#endif