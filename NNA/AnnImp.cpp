#ifndef ANN_IMPLEMENTATION
#define ANN_IMPLEMENTATION
#include "AnnDef.h"
#include <cmath>
#include <random>
#include <iostream>
double NNA::f(Eigen::RowVector2d x) {
	return (100 * std::pow(x[1] - x[0] * x[0], 2.0) + std::pow(1 - x[0], 2.0));
}
void NNA::initialization(Eigen::RowVectorXd &X_LB, Eigen::RowVectorXd &X_UB, Eigen::MatrixXd& x_pattern, Eigen::RowVectorXd& cost) {
	for (size_t i = 0; i < NPOP; i++)
	{
		X_LB [i] = LB;
		X_UB[i] = UB;
	}
	double beta = 1;
	for (size_t i = 0; i < NPOP; i++)
	{
		for (size_t j = 0; j < NVARS; j++) {
			x_pattern(i, j) = 0.0;
		}
	}
	for (size_t i = 0; i < NPOP; i++)
	{
		cost[i] = 0.0;
	}
	//generate randomly initial population
	std::default_random_engine  generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	for (size_t i = 0; i < NPOP; i++)
	{
		for (size_t j = 0; j < NVARS; j++)
		{
			x_pattern(i, j) = LB + ((UB - LB) * distribution(generator));
		}
		cost[i] = f(x_pattern);
	}
}
void NNA::createRandom(Eigen::Matrix<double,1,100>  ww, Eigen::Matrix<double, 100, 100>w) {
	ww = ww.setOnes();
	for (size_t i = 0; i <NPOP; i++)
	{
		ww[i] = ww[i] * 0.5;
	}
	for (size_t i = 0; i < NPOP; i++)
	{
		for (size_t j = 0; j < NPOP; j++) {
			if (i == j) {
				w(i, j) = 1.0;
			}
			else
			{
				w(i, j) = 0.0;
			}
		}
	}
	Eigen::Matrix<double, 1, 99> t; //size t= npop -1
	std::default_random_engine  generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	unsigned int terminat = 0;
	while (terminat!=NPOP)
	{
		for (size_t i = 0; i < NPOP - 1; i++)
		{
			t[i] = distribution(generator) * 0.5;
		}
		std::cout << "phase 1" << std::endl;
		for (size_t i = 0; i < NPOP-1; i++)
		{
			t[i] = t[i] / t.sum() * 0.5;//right division 
		}
		std::cout << "phase 2" << std::endl;

		for (size_t i = 0; i < NPOP; i++)
		{
			for (size_t j = 0; j < NPOP; j++) {
				w(i, j) = 0;
			}
		}
		std::cout << "phase 3" << std::endl;

		++terminat;
	}
	std::cout << "Sum \t: " << ww.sum() << std:: endl;
	std::cout << "Sum t \t: " << t.sum() << std::endl;
	std::cout << "Sum w \t: " << w.sum() << std::endl;

}

#endif 
