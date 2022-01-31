#include <iostream>
#include <Eigen/dense>
#include <random> 
#include <algorithm>
#include "AnnDef.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <fstream>
int n_rotate = 0, n_wrotate = 0;
double target = 0.0, beta = 1.0;
constexpr int FLOAT_MIN = 0, FLOAT_MAX=1;
std::random_device rd; //nondeterminioctic random distribution
std::default_random_engine eng(rd());
std::uniform_real_distribution<float> 
distr(FLOAT_MIN, FLOAT_MAX), distr1(L1, U1), distr2(L2, U2), distr3(L3, U3), distr4(L4, U4);//random number between 1 and 0
std::uniform_int_distribution<int> newDistr(0, nvars - 1), newDistrWeight(0, NPOP - 1);//table 1
Eigen::Matrix <double, 1, 2> value_index1, value_index2;//target value-index

int main()
{

	Eigen::Matrix<double, NPOP, 1>  cost;
	Eigen::Matrix<double, NPOP, nvars> x_pattern, x_new;
	Eigen::Matrix<double, 1, NPOP> ww;
	Eigen::Matrix<double, NPOP, NPOP> w;
	Eigen::Matrix<double, 1, NPOP - 1> t;
	Eigen::Matrix <double, 1, 2> value_index;
	Eigen::Matrix<double, 1, nvars> x_target;
	Eigen::Matrix<double, NPOP, 1> w_target;
	//call NNA
	for (size_t i = 0; i < 20; i++)
	{
		NNA::initialization(ww, w, x_pattern, cost);
		NNA::CreateInitialPopulation(x_pattern, cost, value_index);
		NNA::generateWeight(w, t);
		NNA::setTarget(x_target, x_pattern, value_index, w_target, w);
		NNA::setTarget(x_target, x_pattern, value_index, w_target, w);
		NNA::Run(x_new, x_pattern, w_target, w, cost, x_target);
	}
	std::cout << "finish" << std::endl;
}
double NNA::f(Eigen::Matrix<double, 1, nvars> x) {
	double PENALTY = std::pow(10, 15.0);
	int caseNum = 25;
	double c[NPOP];
	double sumConstrains = 0.0;
	double temp = 0;
	double p = 1.0;
	int n = 2;
	double sum = 0.0;

	//beam design 
	double P = 6000, L = 14, E = 30e+6, G = 12e+6,
	t_max = 13600, s_max = 30000, d_max = 0.25,M,R,J,P_c,t1,t,s,d,t2;
	//note max transform by pre mul-multiplying the objective function by -1
	switch (caseNum)
	{
	case 0:
		//griewank function 
		for (size_t i = 0; i < nvars; i++)
		{
			temp += (std::pow(x(0, i), 2.0) / 4000.0);
		}
		for (size_t i = 1; i < nvars + 1; i++)
		{
			p *= std::cos(x(0, i - 1) / std::sqrt(i));
		}
		return (1.0 + temp - p);
		break;
	case 1:
		//sphere
		temp = 0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp += std::pow(x(0, i), 2.0);
		}
		return temp;
		break;
	case 2:
		//cube  function 
		temp = 0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp += -std::pow(x(0, i), 3.0);
		}
		return temp;
		break;
	case 3:
		//Michalewicz function 
		temp = 0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp += std::sin(x(0, i)) * std::pow((std::sin(i + 1 * std::pow(x(0, i), 2.0) / M_PI)), 2.0 * 10.0);
		}
		return -temp;
		break;
	case 4:
		//sum square
		temp = 0.0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp = temp + (i + 1) * std::pow(x(0, i), 2.0);
		}
		return temp;
		break;
	case 5:
		//rosenbrock function
		temp = 0.0;
		for (size_t i = 1; i < nvars; i++)
		{
			temp = temp + 100.0 * (std::pow(x(0, i) - std::pow(x(0, i - 1), 2.0), 2.0)) + std::pow(1.0 - x(0, i - 1), 2.0);
		}
		return temp;
		break;
	case 6:
		//rastrign function
		temp = 0.0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp = temp + (std::pow(x(0, i), 2.0) - 10.0 * std::cos(2.0 * M_PI * x(0, i)));
		}
		return 10.0 * nvars + temp;
		break;
	case 7:
		//hump function
		temp = 1.0316285 + 4 * std::pow(x(0, 0), 2.0) -
			2.1 * std::pow(x(0, 0), 4.0) + std::pow(x(0, 0), 6.0) / 3 + x(0, 0) * x(0, 1) 
			- 4 * std::pow(x(0, 1), 2.0) + 4 * std::pow(x(0, 1), 4.0);
		return temp;
		break;
	case 8:
		//SCHWEFEL function
		temp = 0.0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp += x(0, i) * std::sin(std::sqrt(std::fabs(x(0, i))));
		}
		return (418.9829 * nvars - temp);
		break;
		//constrained optimization 
	case 9:
		//page 414. stewart calculus : 1.86903e-06 -1.86778e-06
		temp = std::pow(x(0, 0), 2.0) + 2.0 * std::pow(x(0, 1), 2.0);
		sumConstrains = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) - 1.0;
		return (temp + PENALTY * std::pow(std::max(sumConstrains, 0.0), 2.0));
		break;
	case 10:
		//page 417. stewart calculus, question no 19 : 0.707822 0.353192
		temp = std::exp(-x(0, 0) * x(0, 1));
		sumConstrains = std::pow(x(0, 0), 2.0) + 4.0 * std::pow(x(0, 1), 2.0) - 1.0;
		return (temp + PENALTY * std::pow(std::max(sumConstrains, 0.0), 2.0));
		break;
	case 11:
		//simple square case : 1.95418 0.0885032, page 425
		temp = (1 / 2.0) * std::pow(x(0, 0) - 2.0, 2.0) + (1 / 2.0) * std::pow(x(0, 1) - 1 / 2.0, 2.0);
		sumConstrains = -std::pow(x(0, 0) + 1.0, -1.0) + x(0, 1) + (1 / 4.0);
		return (temp + PENALTY * std::pow(std::max(sumConstrains, 0.0), 2.0));
		break;
	case 12:
		//page 475, jorge nocedal  : 1.38838 1.69419
		temp = std::pow(x(0, 0) - 1.0, 2.0) + std::pow(x(0, 1) - 2.5, 2.0);
		c[0] = -x(0, 0) + 2.0 * x(0, 1) - 2.0;
		c[1] = x(0, 0) + 2.0 * x(0, 1) - 6.0;
		c[2] = x(0, 0) - 2.0 * x(0, 1) - 2.0;
		c[3] = 0.0;
		for (size_t i = 0; i < 3; i++)
		{
			c[3] += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + c[3]);
		break;
	case 13:
		// jorge nocedal page 322 : 0.999732 0.000265118, 
		temp = std::pow(x(0, 0) - 3 / 2.0, 2.0) + std::pow(x(0, 1) - 1 / 2.0, 4.0);
		c[0] = -1.0 + x(0, 0) + x(0, 1);
		c[1] = -1.0 + x(0, 0) - x(0, 1);
		c[2] = -1.0 - x(0, 0) + x(0, 1);
		c[3] = -1.0 - x(0, 0) - x(0, 1);
		c[4] = 0.0;
		for (size_t i = 0; i < 5; i++)
		{
			c[4] += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + c[4]);
		break;
	case 14:
		//simple square case : 1.95418 0.0885032, page 425
		temp = (2.0 * std::sqrt(2.0)*x(0,0) + x(0, 1)) * 100.0;
		c[0]=(std::sqrt(2.0) * x(0, 0) + x(0, 1)) / (std::sqrt(2.0) * std::pow(x(0, 0), 2.0) + 2.0 * x(0, 0) * x(0, 1)) * 2.0 - 2.0;
		c[1] = (x(0, 1)) / (std::sqrt(2.0) * std::pow(x(0, 0), 2.0) + 2.0 * x(0, 0) * x(0, 1)) * 2.0 - 2.0;
		c[2] = (1.0) / (x(0, 0) + std::sqrt(2.0)* x(0, 1)) * 2.0 - 2.0;
		return (temp + PENALTY * std::pow(std::max(c[0], 0.0), 2.0) + PENALTY * std::pow(std::max(c[1], 0.0), 2.0) 
			+ PENALTY * std::pow(std::max(c[2], 0.0), 2.0));
		break;
	case 15://max paper zahara  
		temp = -std::pow(std::sin(2.0 * M_PI * x(0, 0)), 3.0) * 
			std::sin(2.0 * M_PI * x(0, 1))/ -(std::pow(x(0,0),3.0) * (x(0,0)+x(0,1)));
		c[0] = std::pow(x(0, 0), 2.0) - x(0, 1) + 1.0;
		c[1] = 1.0 - x(0, 0) + std::pow(x(0, 1) - 4.0, 2.0);
		return (temp + PENALTY * std::pow(std::max(c[0], 0.0), 2.0) + PENALTY * std::pow(std::max(c[1], 0.0), 2.0));
		break;
	case 16:
		return (-x(0, 0) - x(0, 1));
		break;
		//nonlinear integer programming 
	case 17:
		return  std::pow(x(0, 0), 4.0) + std::pow(x(0, 1), 4.0) + 
			16.0 * (x(0, 0) * x(0, 1) + std::pow(4.0 + x(0, 1), 2.0));
		break;
	case 18: 
		temp= 10.0 / x(0, 0) + 10.0 / x(0, 1) + 20.0 / x(0, 2) + 30.0 / x(0, 3);
		c[0] = x(0, 0) + x(0, 1) + x(0, 2) + x(0, 3) - 10.0;
		return (temp + PENALTY * std::pow(std::max(c[0], 0.0), 2.0));
		break;
		//mixed nollinear integer programming 
	case 19:
		temp = -0.7 * x(0, 2) + 5.0 * std::pow((x(0, 0) - 0.5), 2.0) + 0.8;
		c[0] = -std::exp(x(0, 0) - 0.2) - x(0, 1);
		c[1] = x(0, 1) + 1.1*x(0, 2) - 1.0;
		c[2] = x(0, 0) - 1.2 * x(0, 2) - 0.2;
		c[3] = 0.2 - x(0, 0);
		c[4] = x(0, 0) - 1.0;
		c[5] = -2.2554 - x(0, 1);
		c[6] = x(0, 1) + 1.0;
		return (temp + PENALTY * std::pow(std::max(c[0], 0.0), 2.0) + PENALTY * std::pow(std::max(c[1], 0.0), 2.0) 
			+ PENALTY * std::pow(std::max(c[2], 0.0), 2.0) + PENALTY * std::pow(std::max(c[3], 0.0), 2.0)
			+ PENALTY * std::pow(std::max(c[4], 0.0), 2.0) + PENALTY * std::pow(std::max(c[5], 0.0), 2.0)
			+ PENALTY * std::pow(std::max(c[6], 0.0), 2.0));
		break;
	case 20:
		temp = -17 * x(0, 0) - 12.0 * x(0, 1);
		c[0] = 10 * x(0, 0) + 7.0 * x(0, 1) - 40.0;
		c[1] = x(0, 0) + x(0, 1) - 5.0;
		return (temp + PENALTY * std::pow(std::max(c[0], 0.0), 2.0) + PENALTY * std::pow(std::max(c[1], 0.0), 2.0));
		break;
		//enginering problem 
	case 21:
		//tension/ compression spring design problem 
		temp = (x(0, 2) + 2.0) * x(0, 1) * std::pow(x(0, 0), 2.0);
		sum = 0.0;
		c[4] = 71785.0 * std::pow(x(0, 0), 4.0);
		c[5] = std::pow(x(0, 1), 3.0) * x(0, 2);
		c[6] = 4.0 * std::pow(x(0, 1), 2.0) - x(0, 0)*x(0,1);
		c[7] = 12566.0 * (x(0, 1) * std::pow(x(0, 0), 3.0) - std::pow(x(0, 0), 4.0));
		c[8] = (1.0 / 5108.0 * std::pow(x(0, 0), 2.0));
		c[9] = x(0, 0) + x(0, 1);
		c[10] = 140.45 * x(0, 0);
		c[11] = std::pow(x(0, 1), 2.0) * x(0, 2);

		c[0] = 1.0 - (c[5] / c[4]);
		c[1] = ( c[6]/c[7]) + c[8] -(1.0);
		c[2] = 1.0 - ( c[10]/c[11]);
		c[3] = (c[9] / 1.5)-1.0;
		

		for (size_t i = 0; i <= 3; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 22:
		//pressure vessel design 
		temp = 0.6224 * x(0, 0) * x(0, 2) * x(0, 3) + 1.7781*std::pow(x(0,2),2.0)*x(0,1)+
			3.1661*std::pow(x(0,0),2.0)*x(0,3)+ 19.84*std::pow(x(0,1),2.0)*x(0,3);
		c[0] = -x(0, 0) + 0.0193 * x(0, 2);
		c[1] = -x(0, 1) + 0.00945 * x(0, 2);
		c[2] = -M_PI * x(0, 3) * std::pow(x(0, 2), 2.0) - (4.0 / 3.0) * M_PI * std::pow(x(0, 2), 3.0) + 1296000;
		c[3] = x(0, 3) - 240.0;

		 
		for (size_t i = 0; i <= 3; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return temp + sum;
		break;
	case 23:
		//three-bar truss design 
		temp = (2.0 * std::sqrt(2.0) * x(0, 0) + x(0, 1)) * 100.0;
		c[0] = ( (std::sqrt(2.0) * x(0, 0) + x(0, 1)) /(std::sqrt(2.0) * std::pow(x(0, 0),2.0) +2.0*x(0,0)*x(0,1)) )*2.0-2.0;
		c[1] = (( x(0, 1)) / (std::sqrt(2.0) * std::pow(x(0, 0), 2.0) + 2.0 * x(0, 0) * x(0, 1))) * 2.0 - 2.0;
		c[2] = (1.0 / (std::sqrt(2.0) * x(0, 1) + x(0, 0))) * 2.0 - 2.0;
		for (size_t i = 0; i <= 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return temp + sum;
 		break;
	case 24:
		//i-beam design problem 
		temp = 5000.0 / (c[4] + c[1]  + c[2] * std::pow(c[3],2.0)	);

		c[1] = x(0, 0) * std::pow(x(0, 3), 3.0) / 6.0;
		c[2] = 2.0 * x(0, 0) * x(0, 3);
		c[3] = (x(0, 1) - x(0, 3) / 2.0);
		c[4] =x(0,2)*  std::pow (	x(0, 1) - 2.0 * x(0, 3),3.0)	/12.0;

		c[0] = 2.0 * x(0, 0) * x(0, 2) + x(0, 2) * (x(0, 1) - 2.0 * x(0, 3));		 
		for (size_t i = 0; i < 1; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 25:
		//weldead beam design 
		temp = 1.10471 * std::pow(x(0, 0), 2.0) * x(0, 1) + 0.04811 * x(0, 2) * x(0, 3) * (14 + x(0, 1));
		M = P * (L + x(0, 1) / 2.0);
		R = std::sqrt(0.25 * (x(0, 1) * x(0, 1) + (x(0, 0) + std::pow(x(0, 2), 2.0))) );
		J = 2 * (std::sqrt(2) * x(0, 0) * x(0, 1) * (x(0, 1) * x(0, 1) / 12 + 0.25 * (x(0, 0) + x(0, 2) * x(0, 2))) );
		P_c = (4.013 * E / (6 * L * L)) * x(0, 2) * x(0, 3) * x(0, 3) * (1 - 0.25 * x(0, 2) * std::sqrt(E / G) / L);
		t1 = P / (std::sqrt(2) * x(0, 0) * x(0, 1)); 
		t2 = M * R / J;
		t = std::sqrt(t1 * t1 + t1 * t2 * x(0,1) / R + t2 * t2);
		s = 6 * P * L / (x(0, 3) * x(0, 2) * x(0, 2));
		d = 4 * P * std::pow(L,3.0) / (E * x(0,3) * std::pow(x(0,2),3.0));
		c[0] = t - t_max;
		c[1] = s - s_max;;
		c[2] = x(0,0) - x(0,3);
		c[3] = 0.10471 * std::pow(x(0,0),2.0) + 0.04811 * x(0,2) * x(0,3) * (14.0 + x(0,1)) - 5.0;
		c[4] = 0.125 - x(0, 0);
		c[5] = d - d_max;
		c[6] = P - P_c;
		for (size_t i = 0; i <= 6; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 26:
		//speed reducer design 
		temp = 0.7854 * x(0, 0) * std::pow(x(0, 1), 2.0) *
			(3.3333*std::pow(x(0,2),2.0)+14.9334*x(0,2)-43.0934)  - 
			1.508*x(0,0) *(std::pow(x(0,5),2.0)+std::pow(x(0,6),2.0))
			+ 7.4777*(std::pow(x(0,5),3.0) + std::pow(x(0, 6), 3.0)) + 
			0.7854*(x(0,3) *std::pow(x(0,5),2.0) + x(0,4)* std::pow(x(0,6),2.0));
		c[0] = (	27.0  / (x(0, 0) * std::pow(x(0, 1), 2.0) * x(0, 2))) -1.0;
		c[1] = (	397.5 / (x(0, 0) * std::pow(x(0, 1), 2.0) * std::pow(x(0, 2),2.0) )) - 1.0;
		c[2] = (	1.93*std::pow(x(0,3),4.0) / (x(0, 1) * std::pow(x(0, 5), 4.0) * x(0, 2))) - 1.0;
		c[3] = (	1.93 * std::pow(x(0, 4), 3.0) / (x(0, 1) * std::pow(x(0, 6), 4.0) * x(0, 2))) - 1.0;
		
		c[4] = 1.0 / 110 * std::pow(x(0, 5), 3.0)  * std::sqrt(	std::pow(745*x(0,3)/x(0,1)*x(0,2),2.0)+16.9*1e6	)-1.0;
		c[5] = 1.0 / 85.0 * std::pow(x(0, 6), 3.0) * std::sqrt(std::pow(745 * x(0, 4) / x(0, 1) * x(0, 2), 2.0) +
		157.5*1e6)-1.0;
		c[6] = (x(0, 1) * x(0, 2) / 40.0) - 1.0;
		c[7] = (x(0, 0) / 12.0 * x(0, 1)) - 1.0;
		c[8] = (5.0 * x(0, 1) / x(0, 0)) - 1.0;
		c[9] = ((1.5 * x(0, 5) + 1.9) / x(0, 3)) - 1.0;
		c[10] = ((1.1 * x(0, 6) + 1.9) / x(0, 4)) - 1.0;
		for (size_t i = 0; i <= 10; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return temp + sum;
		break;
	//multiobjective optimization 
	case 27:
		c[0] = std::pow(x(0, 0),2.0);
		c[1] = std::pow(x(0, 0) - 2.0, 2.0);
		for (size_t i = 0; i < 2; i++)
		{
			sum += 0.5 * c[i];
		}
		return sum;
		break;
	default:
		break;
	}
	return 0;
}
void NNA::initialization(Eigen::Matrix<double, 1, NPOP>& ww, Eigen::Matrix<double, NPOP, NPOP>& w,
	Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& cost) {
	ww = ww.setOnes() * 0.5;
	w = 0.5 * w.setIdentity();
	x_pattern = x_pattern.setZero();
	cost = cost.setZero();
}
void NNA::CreateInitialPopulation(Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& cost,
	Eigen::Matrix <double, 1, 2>& value_index) {
	for (size_t i = 0; i < NPOP; i++)
	{
		for (size_t j = 0; j < nvars; j++)
		{
			//modify
			//tension 
			/*
			f ((j == 0))
			{
				//x_pattern(i, j) = L1 + (U1 - L1) * distr(eng);
				x_pattern(i, j) = distr1(eng);
			}
			else if (j == 1) {
				//x_pattern(i, j) = L2 + (U2 - L2) * distr(eng);
				x_pattern(i, j) = distr2(eng);
			}
			else {
				//x_pattern(i, j) = L3 + (U3 - L3) * distr(eng);
				x_pattern(i, j) = distr3(eng);

			}
			*/

			/* pressure vessel 
			if ((j == 0) || (j == 1))
			{
				//x_pattern(i, j) = L1 + (U1 - L1) * distr(eng);
				x_pattern(i, j) = distr1(eng);
			}
			else {
				//x_pattern(i, j) = L3 + (U3 - L3) * distr(eng);
				x_pattern(i, j) = distr3(eng);

			}*/

			// i-beam design 
			/*if ((j == 2) || (j == 3))
			{
				//x_pattern(i, j) = L1 + (U1 - L1) * distr(eng);
				x_pattern(i, j) = distr3(eng);
			}
			else if (j==1)
			{
				x_pattern(i, j) = distr2(eng);
			}
			else {
				//x_pattern(i, j) = L3 + (U3 - L3) * distr(eng);
				x_pattern(i, j) = distr1(eng);
			}
			*/
			if ((j == 0) || (j == 3))
			{
				//x_pattern(i, j) = L1 + (U1 - L1) * distr(eng);
				x_pattern(i, j) = distr1(eng);
			}
			else {
				x_pattern(i, j) = distr2(eng);

			}
			//stop
			//x_pattern(i, j) = LB + (UP - LB) * distr(eng);
		}
		cost(i, 0) = NNA::f(x_pattern.row(i));
	}
	for (size_t i = 0; i < cost.rows(); i++)
	{
		if (cost(i, 0) == cost.minCoeff()) { //min
			value_index(0, 0) = cost(i, 0);
			value_index(0, 1) = i;
			break;
		}
	}
}
void NNA::generateWeight(Eigen::Matrix<double, NPOP, NPOP>& w, Eigen::Matrix<double, 1, NPOP - 1>& t) {
	for (int var = 0; var < w.cols(); var++)
	{
		double total = 0.0;

		for (int j = 0; j < t.cols(); j++)
		{
			t(0, j) = distr(eng) * 0.5;
			total += t(0, j);

		}
		for (int k = 0; k < t.cols(); k++)
		{
			t(0, k) = (t(0, k) / total) * 0.5;
		}
		int ink = 0;
		for (size_t j = 0; j < w.cols(); j++)
		{
			if (w(j, var) == 0.0) {
				w(j, var) = t(ink);
				++ink;
			}
		}
	}
}
void NNA::setTarget(Eigen::Matrix<double, 1, nvars>& x_target, Eigen::Matrix<double, NPOP, nvars>& x_pattern,
	Eigen::Matrix <double, 1, 2>& value_index, Eigen::Matrix<double, NPOP, 1>& w_target,
	Eigen::Matrix<double, NPOP, NPOP>& w) {
	x_target = x_pattern.row(int(value_index(0, 1)));
	target = value_index(0, 0);
	w_target = w.col(int(value_index(0, 1)));
}
void NNA::Run(Eigen::Matrix<double, NPOP, nvars>& x_new, Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& w_target,
	Eigen::Matrix<double, NPOP, NPOP>& w, Eigen::Matrix<double, NPOP, 1>& cost, Eigen::Matrix<double, 1, nvars>& x_target) {
	std::ofstream writeF;
	auto Maximum = 151;
	int begin = 1;
	int auxilaryVar = 0;
	writeF.open("weldeadBeam.csv", std::ios::app);
	while (begin != Maximum)
	{
		//step 6: generate new pattern and update solution... 
		x_new = w * x_pattern;
		x_pattern = x_new + x_pattern;
		//return new_solution.stop 
		//step 7: update the weight matrix..
		for (size_t i = 0; i < NPOP; i++)
		{
			for (size_t j = 0; j < NPOP; j++)
			{
				w(j, i) = std::fabs(w(j, i) + (2.0 * distr(eng)) * (w_target(j, 0) - w(j, i)));
			}
		}
		//consider the constrains for update the weight 
		for (size_t i = 0; i < NPOP; i++)
		{
			w.col(i) = (1 / w.col(i).sum()) * w.col(i);
		}
		//return new weight.... stop

		//step 8: check the bias condition, then perform the bias operator... 
		//procedure in table 1
		for (size_t i = 0; i < NPOP; i++)
		{
			if (distr(eng) < beta) {
				n_rotate = std::ceil(nvars * beta);

				for (size_t l = 0; l < n_rotate; l++)
				{
					//modify here
					auxilaryVar = newDistr(eng);
					//tension
					/*
					if (2 == auxilaryVar)
					{
						//x_pattern(i, auxilaryVar) = L3 + (U3 - L3) * distr(eng);
						x_pattern(i, auxilaryVar) = distr3(eng);
					}
					else if (1== auxilaryVar)
					{
						//x_pattern(i, auxilaryVar) = L2 + (U2 - L2) * distr(eng);
						x_pattern(i, auxilaryVar) = distr2(eng);
					}
					else {
						//x_pattern(i, auxilaryVar) = L1 + (U1 - L1) * distr(eng);
						x_pattern(i, auxilaryVar) = distr1(eng);
					}
					*/ 
					//pressure Vessel
					/*
					if ((auxilaryVar==1)|| (auxilaryVar==2))
					{
						//x_pattern(i, auxilaryVar) = L3 + (U3 - L3) * distr(eng);
						x_pattern(i, auxilaryVar) = distr1(eng);
					}
					else {
						x_pattern(i, auxilaryVar) = distr2(eng); 
					}
					*/
					//i-Beam 
					/*
					if ((auxilaryVar == 2) || (auxilaryVar == 3))
					{
						//x_pattern(i, auxilaryVar) = L3 + (U3 - L3) * distr(eng);
						x_pattern(i, auxilaryVar) = distr3(eng);
					}
					else if (auxilaryVar==1)
					{
						x_pattern(i, auxilaryVar) = distr2(eng);
					}
					else {
						x_pattern(i, auxilaryVar) = distr1(eng);
					}
					*/
					if ((auxilaryVar == 0) || (auxilaryVar == 3))
					{
						x_pattern(i, auxilaryVar) = distr1(eng);
					}
					else {
						x_pattern(i, auxilaryVar) = distr2(eng);
					}
					// stop
				}
				n_wrotate = std::ceil(NPOP * beta);
				for (size_t m = 0; m < n_wrotate; m++)
				{
					w(m, newDistrWeight(eng)) = distr(eng);
				}
				//apply constrains for weoght matrix
				for (size_t n = 0; n < NPOP; n++)
				{
					w.col(n) = (1 / w.col(n).sum()) * w.col(n);//reformulate for element wise division
				}
			}
			//return new weight matrix and x_pattern, with a restriction that rand<beta. stop 

			//step 9, a mistake step in paper, written step 8 : new step : 9
			//apply transfer function opearator, formulate placed in table 2
			else {
				for (size_t p = 0; p < nvars; p++)
				{
					x_pattern(i, p) = x_pattern(i, p) + (2.0 * distr(eng)) * (x_target(0, p) - x_pattern(i, p));
					//equation 13, page 751
					//add auxilary makro for mixed integer problem. test case problem 19
					/*if (p == 0 || p == 1)//position of the integer problem. for 
					{
						x_pattern(i, p) = x_pattern(i, p) + (2.0 * distr(eng)) * (std::round(x_target(0, p)) - x_pattern(i, p));
					}*/
				}
			}
			//return x_pattern, which is satisfy the restriction. stop 
		}
		//step 10: calculate the objective function 
		//before calculate, consider the lower and upper bound
		for (size_t x_max = 0; x_max < NPOP; x_max++)
		{
			for (size_t y_max = 0; y_max < nvars; y_max++)
			{
				//modify here
				//tension ...
				/*
				if (y_max == 0)
				{
					if ( (x_pattern(x_max, y_max) < L1 ) || (x_pattern(x_max, y_max) > U1) )
					{
						x_pattern(x_max, y_max) = distr1(eng);
					}
				}
				else if (y_max == 1) {
					if ((x_pattern(x_max, y_max) < L2) || (x_pattern(x_max, y_max) > U2))
					{

						x_pattern(x_max, y_max) = distr2(eng);
					}
				}
				else {
					if ((x_pattern(x_max, y_max) < L3) || (x_pattern(x_max, y_max) > U3))
					{
						x_pattern(x_max, y_max) = distr3(eng);
					}
				}
				*/
				//pressure Vessel design 
				/*
				if ((y_max == 0) || (y_max==1))
				{
					if ((x_pattern(x_max, y_max) < L1) || (x_pattern(x_max, y_max) > U1))
					{
						x_pattern(x_max, y_max) = distr1(eng);
					}
				}
				else
				{
					if ((x_pattern(x_max, y_max) < L3) || (x_pattern(x_max, y_max) > U3))
					{
						x_pattern(x_max, y_max) = distr3(eng);
					}
				}
				*/
				/*
				* //ibeam design 
				if ((y_max == 2) || (y_max == 3))
				{
					if ((x_pattern(x_max, y_max) < L3) || (x_pattern(x_max, y_max) > U3))
					{
						x_pattern(x_max, y_max) = distr3(eng);
					}
				}
				else if (y_max==1)
				{
					if ((x_pattern(x_max, y_max) < L2) || (x_pattern(x_max, y_max) > U2))
					{
						x_pattern(x_max, y_max) = distr2(eng);
					}
				}
				else
				{
					if ((x_pattern(x_max, y_max) < L1) || (x_pattern(x_max, y_max) > U1))
					{
						x_pattern(x_max, y_max) = distr1(eng);
					}
				}
				*/
				if ((y_max == 0) || (y_max == 3))
				{
					if ((x_pattern(x_max, y_max) < L1) || (x_pattern(x_max, y_max) > U1))
					{
						x_pattern(x_max, y_max) = distr1(eng);
					}
				}
				else 
				{
					if ((x_pattern(x_max, y_max) < L2) || (x_pattern(x_max, y_max) > U2))
					{
						x_pattern(x_max, y_max) = distr2(eng);
					}
				}
			}
		}

		
		//when calculate the cost, we already assume that all value not outside the lower and upper bound
		for (size_t q = 0; q < NPOP; q++)
		{
			cost(q, 0) = NNA::f(x_pattern.row(q));
		}
		//return all cost. stop 

		//step 11: update the target solution 
		//find minimum cost
		for (size_t r = 0; r < NPOP; r++)
		{
			if (cost(r, 0) == cost.minCoeff()) {//min coeff- transform the proble to max
				value_index1(0, 0) = cost(r, 0);
				value_index1(0, 1) = r;
				break;
			}
		}

		//find max cost
		for (size_t s = 0; s < NPOP; s++)
		{
			if (cost(s, 0) == cost.maxCoeff()) {//max coeff.- transform the proble to min
				value_index2(0, 0) = cost(s, 0);
				value_index2(0, 1) = s;
				break;
			}
		}
		if (value_index1(0, 0) < target)
		{
			target = value_index1(0, 0);
			x_target = x_pattern.row(int(value_index1(0, 1)));
			w_target = w.col(int(value_index1(0, 1)));
		}
		else {
			x_pattern.row(int(value_index2(0, 1))) = x_target;
			w.col(int(value_index2(0, 1))) = w_target;
		}
		//step 12: update beta, equation 11
		beta *= 0.99;
		if (beta < 0.01)
		{
			beta = 0.05;
		}
		//std::cout << target << std::endl;

		//return new beta. stop 
		writeF << begin << "\t;";
		for (size_t i = 0; i < nvars; i++)
		{
			writeF << x_target(0,i) << "\t;";

		}
		writeF << target << std::endl;   
		
		//std::cout << x_target << "\t" << target << std::endl;	
		++begin;
	}
	writeF << std::endl;
	writeF.close();
	//std::cout << x_target << std::endl;

	//std::cout << target << std::endl;
	//std::cout << "xtarge \t: " <<x_target << std::endl;
}
