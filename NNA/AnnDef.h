#pragma once
#ifndef ANN_HEADER
#define ANN_HEADER
#include <Eigen/dense>
#define LB 12.0
#define UB 15.0
#define NVARS 2
#define NPOP 100
#define MAX_ITER 1000
namespace NNA {
	static double f(Eigen::RowVector2d);
	static void initialization(Eigen::RowVectorXd &,Eigen::RowVectorXd &, Eigen::MatrixXd&,Eigen::RowVectorXd&);
	static void createRandom(Eigen::Matrix<double, 1, 100>  , Eigen::Matrix<double, 100, 100>);
	static void constrainsSummation();
	static void bestSolution();
	static void transferOperator();
	static void biasReduction();
	static void biasForWeight();
	static void updateWeight();
	static Eigen::RowVector2d mainNNA();
}
#endif 

