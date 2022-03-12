#include <iostream>
#include <Eigen/dense>
#include <random> 
#include <algorithm>
#include "AnnDef.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <chrono>
#define KODE 29
#define FKODE 29
#define maxSamp  1500//case monte carlo
int n_rotate = 0, n_wrotate = 0;
double target = 0.0, beta = 1.0;
constexpr int FLOAT_MIN = 0, FLOAT_MAX=1;
std::random_device rd; //nondeterminioctic random distribution
std::default_random_engine eng(rd());
std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);
//std::uniform_real_distribution<double> boundMonteA(0,M_PI/2.0), boundMonteB(0,1.0);
std::uniform_real_distribution<double> boundMonteA(0,3), boundMonteB(1,2),monteTripA(0,1);

//std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX), distr1(L1, U1), distr2(L2, U2), distr3(L3, U3), distr4(L4, U4);//random number between 1 and 0
std::uniform_real_distribution<double>
d1(LoSpeedRed[0], UpSpeedRed[0]),
d2(LoSpeedRed[1], UpSpeedRed[1]),
d3(LoSpeedRed[2], UpSpeedRed[2]),
d4(LoSpeedRed[3], UpSpeedRed[3]),
d5(LoSpeedRed[4], UpSpeedRed[4]),
d6(LoSpeedRed[5], UpSpeedRed[5]),
d7(LoSpeedRed[6], UpSpeedRed[6]), cb(cantBeam[0], cantBeam[1]),
cur1(currugatedBulkhead[0], currugatedBulkhead[1]), cur2(currugatedBulkhead[2], currugatedBulkhead[3]),
tC(tabColumn[0], tabColumn[1]), ge(gear[0], gear[1]),
md1(LoMulD[0], UpMulD[0]), md2(LoMulD[1], UpMulD[1]), md3(LoMulD[2], UpMulD[2]), md4(LoMulD[3], UpMulD[3]),
md5(LoMulD[4], UpMulD[4]), 
tn1(tension[0], tension[1]), tn2(tension[2], tension[3]), tn3(tension[4], tension[5]),
pr1(pressureVes[0] , pressureVes[1]), pr2(pressureVes[2], pressureVes[3]), pr3(pressureVes[4], pressureVes[5]), 
pr4(pressureVes[6], pressureVes[7]), 
ib1(ibeam[0] , ibeam[1]), ib2(ibeam[2], ibeam[3]), ib3(ibeam[4], ibeam[5]), ib4(ibeam[6], ibeam[7]),
wb1(weldeadBeam[0] , weldeadBeam[1]) , wb2(weldeadBeam[2], weldeadBeam[3]);
std::uniform_int_distribution<int> newDistr(0, nvars - 1), newDistrWeight(0, NPOP - 1) , 
muldiskF(0,40) , muldiskT(0,3);//table 1
//defenition for default problem without specific bound
std::uniform_real_distribution<double> DF(def[0], def[1]),
CS1(cs1421[0] , cs1421[1]), CS2(cs1421[2], cs1421[3]), 
DS1(cs15[0], cs15[1]), DS2(cs15[2], cs15[3]), DS3(cs15[4], cs15[5]),
ES1(cs16[0], cs16[1]), ES2(cs16[2], cs16[3]), ES3(cs16[4], cs16[5]),
FS1(cs19[0], cs19[1]), FS2(cs19[2], cs19[3]), FS3(cs19[4], cs19[5]), 
GS1(cs21[0], cs21[1]), GS2(cs21[2], cs21[3]);
std::uniform_real_distribution<double>CS29_32(cs29_32[0], cs29_32[1]), DS29_32(cs29_32[2], cs29_32[3]);

Eigen::Matrix <double, 1, 2> value_index1, value_index2;//target value-index
double fx(Eigen::Matrix<double,1,2> x);

int main()
{

	Eigen::Matrix<double, NPOP, 1>  cost;
	Eigen::Matrix<double, NPOP, nvars> x_pattern, x_new;
	Eigen::Matrix<double, 1, NPOP> ww;
	Eigen::MatrixXd w(NPOP,NPOP);
	Eigen::Matrix<double, 1, NPOP - 1> t;
	Eigen::Matrix <double, 1, 2> value_index;
	Eigen::Matrix<double, 1, nvars> x_target;
	Eigen::Matrix<double, NPOP, 1> w_target;
	//case multiple disk 
	for (size_t i = 1; i < 41; i++)
	{
		bF[i] = bF[0] + (i * 10.0);
	}
	//call NNA
	for (size_t i = 0; i < 1; i++)
	{
		NNA::initialization(ww, w, x_pattern, cost);
		NNA::CreateInitialPopulation(x_pattern, cost, value_index);
		NNA::generateWeight(w, t);
		NNA::setTarget(x_target, x_pattern, value_index, w_target, w);
		NNA::setTarget(x_target, x_pattern, value_index, w_target, w);
		NNA::Run(x_new, x_pattern, w_target, w, cost, x_target);
	}

	//monte carlo test 
	/*
	Eigen::MatrixXd  xy(maxSamp,2),z(maxSamp,1), fValue(maxSamp,1);
	
	double area = M_PI / 2.0;
	
	for (size_t i = 0; i < maxSamp; i++)
	{
		z(i, 0) = boundMonteA(eng);
	}
	for (size_t i = 0; i < maxSamp; i++)
	{
		xy(i, 0) = z(i,0);
	}
	for (size_t i = 0; i < maxSamp; i++)
	{
		z(i, 0) = boundMonteB(eng);
	}
	for (size_t i = 0; i < maxSamp; i++)
	{
		xy(i, 1) = z(i, 0);
	}
	for (size_t i = 0; i < maxSamp; i++)
	{
		fValue(i, 0) = fx(xy.row(i));
	}
	
	std::cout << "estimate monte carlo for double integral \t: " << M_PI / 2.0 * fValue.mean(); 
	//std::cout << "estimate monte carlo for double integral \t: " << 3.0 * fValue.mean();
	*/
	//estimate triple integral 
	/*
	double startT, finishT;
	startT = clock();
	double count = 0;
	int maxIter =maxSamp;
	Eigen::MatrixXd sample(maxIter, 3);
	Eigen::MatrixXd temp(maxIter, 1);
	for (size_t i = 0; i < maxSamp; i++)
	{
		temp(i, 0) = std::pow(monteTripA(eng),2.0);
	}
	sample.col(0) = temp.col(0);
	for (size_t i = 0; i < maxSamp; i++)
	{
		temp(i, 0) = std::pow(monteTripA(eng), 2.0);
	}
	sample.col(1) = temp.col(0);
	for (size_t i = 0; i < maxSamp; i++)
	{
		temp(i, 0) = std::pow(monteTripA(eng), 2.0);
	}
	sample.col(2) = temp.col(0);

	for (size_t i = 0; i < maxSamp; i++)
	{
		if (sample.row(i).sum()<1)
		{
			++count;
		}
	}
	finishT = clock();
	std::cout << "\n estimate monte carlo for triple integral \t: " << count / maxSamp;
	std::cout << "\nTime\t: " << (finishT - startT)/CLOCKS_PER_SEC << std::endl;
	

/*	Eigen::Matrix < double, 1, nvars> x;
	double c[10];
	double temp;
	x << 2.32952, 3.17849;
	temp = -x(0, 0) - x(0, 1);
	c[0] = x(0, 1) - 2.0 * std::pow(x(0, 0), 4.0)
		+ 8.0 * std::pow(x(0, 0), 3.0) - 8.0 * std::pow(x(0, 0), 2.0) - 2.0;
	c[1] = x(0, 1) - 4.0 * std::pow(x(0, 0), 4.0)
		+ 32 * std::pow(x(0, 0), 3.0) - 88.0 * std::pow(x(0, 0), 2.0) + 96.0 * x(0, 0) - 36.0;
	std::cout << temp << std::endl;
	for (size_t i = 0; i < 2; i++)
	{
		std::cout << c[i] << std::endl;
	}
	
	* /
	/*
	x << 0.724348, 0.398977;
	temp = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) - std::cos(17.0 * x(0, 0)) - std::cos(17.0 * x(0, 1)) + 3.0;
	c[0] = std::pow(x(0, 0) - 2.0, 2.0) + std::pow(x(0, 1), 2.0) - std::pow(1.6, 2.0);
	c[1] = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1) - 3.0, 2.0) - std::pow(2.7, 2.0);
	std::cout << temp << std::endl;
	for (size_t i = 0; i < 2; i++)
	{
		std::cout << c[i] << std::endl;
	} */

	/*x << 100.017, 1732.08, 6245.52, 45.0394, 250.345, 106.99, 193.238, 350.334;
	temp = x(0, 0) + x(0, 1) + x(0, 2);
	c[0] = -1.0 + 0.0025 * (x(0, 3) + x(0, 5));
	c[1] = -1.0 + 0.0025 * (x(0, 4) + x(0, 6) - x(0, 3));
	c[2] = -1 + 0.01 * (x(0, 7) - x(0, 4));
	c[3] = -x(0, 0) * x(0, 5) + 833.33252 * x(0, 3) + 100 * x(0, 0) - 83333.333;
	c[4] = -x(0, 1) * x(0, 6) + 1250 * x(0, 4) + x(0, 1) * x(0, 3) - 1250 * x(0, 3);
	c[5] = -x(0, 2) * x(0, 7) + 1250000 + x(0, 2) * x(0, 4) - 2500 * x(0, 4);
	std::cout << temp << std::endl;
	for (size_t i = 0; i < 6; i++)
	{
		std::cout << c[i] << std::endl;
	} */
	/* 
	x << 0.689317 ,0.844659;
	temp = std::pow(x(0, 0) - 2.0, 2.0) + std::pow(x(0, 1) - 1.0, 2.0);
	c[0] = x(0, 0) - 2.0 * x(0, 1) + 1.0;//equality 
	c[1] = (std::pow(x(0, 0), 2.0) / 4.0) + std::pow(x(0, 1), 2.0) - 1.0;
	std::cout << temp << std::endl;
	for (size_t i = 0; i < 2; i++)
	{
		std::cout << c[i] << std::endl;
	}*/

	/*
	x << 1.98115 , 2.59376 , 8.91235 , 5.65081, 0.923705 , 1.21431  ,1.08751 , 9.70767  ,7.61723 , 7.71233;
	temp = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) + x(0, 0) * x(0, 1) - 14.0 * x(0, 0) - 16.0 * x(0, 1) +
		std::pow(x(0, 2) - 10.0, 2.0) + 4.0 * std::pow(x(0, 3) - 5.0, 2.0) + std::pow(x(0, 4) - 3.0, 2.0) +
		2.0 * std::pow(x(0, 5) - 1.0, 2.0) + 5.0 * std::pow(x(0, 6), 2.0) + 7.0 * std::pow(x(0, 7) - 11.0, 2.0) +
		2.0 * std::pow(x(0, 8) - 10.0, 2.0) + std::pow(x(0, 9) - 7.0, 2.0) + 45.0;
	c[0] = -105.0 + 4.0 * x(0, 0) + 5.0 * x(0, 1) - 3.0 * x(0, 6) + 9.0 * x(0, 7);
	c[1] = 10 * x(0, 0) - 8.0 * x(0, 1) - 17.0 * x(0, 6) + 2.0 * x(0, 7);
	c[2] = -8.0 * x(0, 0) + 2.0 * x(0, 1) + 5.0 * x(0, 8) - 2.0 * x(0, 9) - 12.0;
	c[3] = 3.0 * std::pow(x(0, 0) - 2.0, 2.0) + 4.0 * std::pow(x(0, 1) - 3.0, 2.0) + 2.0 * std::pow(x(0, 2), 2.0) - 7.0 * x(0, 3) - 120.0;
	c[4] = 5.0 * std::pow(x(0, 0), 2.0) + 8.0 * x(0, 1) + std::pow(x(0, 2) - 6.0, 2.0) - 2.0 * x(0, 3) - 40.0;
	c[5] = std::pow(x(0, 0), 2.0) + 2.0 * std::pow(x(0, 1) - 2.0, 2.0) - 2.0 * x(0, 0) * x(0, 1) + 14.0 * x(0, 4) - 6.0 * x(0, 5);
	c[6] = 0.5 * std::pow(x(0, 0) - 8.0, 2.0) + 2.0 * std::pow(x(0, 1) - 4.0, 2.0) + 3.0 * std::pow(x(0, 4), 2.0) - x(0, 5) - 30.0;
	c[7] = -3.0 * x(0, 0) + 6.0 * x(0, 1) + 12.0 * std::pow(x(0, 8) - 8.0, 2.0) - 7.0 * x(0, 9);

	std::cout << temp << std::endl;
	for (size_t i = 0; i < 8; i++)
	{
		std::cout << c[i] << std::endl;
	}
	*/
	/*
	x << 676.219, 1030, 0.121521, -0.394958;
	temp = 3.0 * x(0, 0) + 0.000001 * std::pow(x(0, 0), 2.0) + 2.0 * x(0, 1)
		+ (0.000002 / 3.0) * std::pow(x(0, 1), 3.0);
	c[0] = -x(0, 3) + x(0, 2) - 0.55; //inequality 
	c[1] = -x(0, 2) + x(0, 3) - 0.55;
	c[2] = 1000.0 * std::sin(-x(0, 2) - 0.25) + 1000.0 * std::sin(-x(0, 3) - 0.25) + 894.8 - x(0, 0);
	c[3] = 1000.0 * std::sin(x(0, 2) - 0.25) + 1000.0 * std::sin(x(0, 2) - x(0, 3) - 0.25) + 894.8 - x(0, 1);
	c[4] = 1000.0 * std::sin(x(0, 3) - 0.25) + 1000.0 * std::sin(x(0, 3) - x(0, 2) - 0.25) + 1294.8;
	x << 853.217, 847.062 ,666.016, 646.713;
	std::cout << temp << std::endl;
	for (size_t i = 0; i < 5; i++)
	{
		std::cout << c[i] << std::endl;
	}
	*/
	/*
	x << 2.24678 ,2.38107;
	temp = std::pow(std::pow(x(0, 0), 2.0) + x(0, 1) - 11.0, 2.0) + std::pow(std::pow(x(0, 1), 2.0) + x(0, 0) - 7.0, 2.0);
	c[0] = std::pow(x(0, 0) - 0.05, 2.0) + std::pow(x(0, 1) - 2.5, 2.0) - 4.84;
	c[1] = -std::pow(x(0, 0), 2.0) - std::pow(x(0, 1) - 2.5, 2.0) + 4.84;
	std::cout << temp << std::endl;
	for (size_t i = 0; i < 2; i++)
	{
		std::cout << c[i] << std::endl;
	}
	*/
	/*
	x << 2.15726 ,  1.96219, - 0.65946,   4.38907 ,- 0.635137   ,1.11285  , 1.47276;
	double temp = std::pow(x(0, 0) - 10.0, 2.0) + 5.0 * std::pow(x(0, 1) - 12.0, 2.0) + std::pow(x(0, 2), 4.0) +
		3.0 * std::pow(x(0, 3) - 11.0, 2.0) + 10.0 * std::pow(x(0, 4), 6.0) + 7.0 * std::pow(x(0, 5), 2.0) +
		std::pow(x(0, 6), 4.0) - 4.0 * x(0, 5) * x(0, 6) - 10.0 * x(0, 5) - 8.0 * x(0, 6);
	c[1] = -127.0 + 2.0 * std::pow(x(0, 0), 2.0) + 3.0 * std::pow(x(0, 1), 4.0) + x(0, 2) + 4.0 * std::pow(x(0, 3), 2.0) +
		5.0 * x(0, 4);
	c[2] = -282.0 + 7.0 * x(0, 0) + 3.0 * x(0, 1) + 10.0 * std::pow(x(0, 2), 2.0) + x(0, 3) - x(0, 4);
	c[3] = -196.0 + 23.0 * x(0, 0) + std::pow(x(0, 1), 2.0) + 6.0 * std::pow(x(0, 5), 2.0) - 8.0 * x(0, 6);
	c[4] = 4.0 * std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) - 3.0 * x(0, 0) * x(0, 1) +
		2.0 * std::pow(x(0, 2), 2.0) + 5.0 * x(0, 5) - 11.0 * x(0, 6);
	std::cout << temp << std::endl;

	for (size_t i = 1; i < 5; i++)
	{
		std::cout << c[i] << std::endl;
	}
	*/

	/* << 78.0221, 35.0661, 31.3669, 44.989, 33.4575;
	double temp = 5.3578547 * std::pow(x(0, 2), 2.0) + 0.8356891 * x(0, 0) * x(0, 4) +
		37.293239 * x(0, 0) - 40729.141;
	c[0] = 85.334407 + 0.0056858 * x(0, 1) * x(0, 4) + 0.0006262 * x(0, 0) * x(0, 3)
		- 0.0022053 * x(0, 2) * x(0, 4) - 92.0;
	c[1] = -85.334407 - 0.0056858 * x(0, 1) * x(0, 4) - 0.0006262 * x(0, 0) * x(0, 3)
		- 0.0022053 * x(0, 2) * x(0, 4) + 92.0;
	c[2] = 80.51249 + 0.0071317 * x(0, 1) * x(0, 4) + 0.0029955 * x(0, 0) * x(0, 1)
		+ 0.0021813 * std::pow(x(0, 2), 2.0) - 110.0;
	c[3] = -80.51249 - 0.0071317 * x(0, 1) * x(0, 4) - 0.0029955 * x(0, 0) * x(0, 1)
		- 0.0021813 * std::pow(x(0, 2), 2.0) + 90.0;
	c[4] = 9.300961 + 0.0047026 * x(0, 2) * x(0, 4) + 0.0012547 * x(0, 0) * x(0, 2)
		+ 0.0019085 * x(0, 2) * x(0, 3) - 25.0;
	c[5] = -9.300961 - 0.0047026 * x(0, 2) * x(0, 4) - 0.0012547 * x(0, 0) * x(0, 2)
		- 0.0019085 * x(0, 2) * x(0, 3) + 20.0;
	std::cout << temp << std::endl;
	for (size_t i = 0; i < 6; i++)
	{
		std::cout << c[i] << std::endl;
	} */


}
double NNA::f(Eigen::Matrix<double, 1, nvars> x) {
	double PENALTY = std::pow(10, 15.0);
	int caseNum = 31;
	double c[NPOP],dl[NPOP];
	double sumConstrains = 0.0;
	double temp = 0;
	double p = 1.0;
	int n = 2;
	double sum = 0.0;
	std::vector<double> mC(14) ;
	//beam design 
	double P = 6000, L = 14, E = 30e+6, G = 12e+6,
	t_max = 13600, s_max = 30000, d_max = 0.25,M,R,J,P_c,t1,t,s,d,t2;
	//note max transform by pre mul-multiplying the objective function by -1
	switch (FKODE)
	{
	case 0:
		//tension/ compression spring design problem 
		temp = (x(0, 2) + 2.0) * x(0, 1) * std::pow(x(0, 0), 2.0);
		sum = 0.0;
		c[4] = 71785.0 * std::pow(x(0, 0), 4.0);
		c[5] = std::pow(x(0, 1), 3.0) * x(0, 2);
		c[6] = 4.0 * std::pow(x(0, 1), 2.0) - x(0, 0) * x(0, 1);
		c[7] = 12566.0 * (x(0, 1) * std::pow(x(0, 0), 3.0) - std::pow(x(0, 0), 4.0));
		c[8] = (1.0 / 5108.0 * std::pow(x(0, 0), 2.0));
		c[9] = x(0, 0) + x(0, 1);
		c[10] = 140.45 * x(0, 0);
		c[11] = std::pow(x(0, 1), 2.0) * x(0, 2);

		c[0] = 1.0 - (c[5] / c[4]);
		c[1] = (c[6] / c[7]) + c[8] - (1.0);
		c[2] = 1.0 - (c[10] / c[11]);
		c[3] = (c[9] / 1.5) - 1.0;


		for (size_t i = 0; i <= 3; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 1:
		//pressure vessel design 
		temp = 0.6224 * x(0, 0) * x(0, 2) * x(0, 3) + 1.7781 * std::pow(x(0, 2), 2.0) * x(0, 1) +
			3.1661 * std::pow(x(0, 0), 2.0) * x(0, 3) + 19.84 * std::pow(x(0, 1), 2.0) * x(0, 3);
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
	
	case 2:
		//i-beam design problem 
		temp = 5000.0 / (c[4] + c[1] + c[2] * std::pow(c[3], 2.0));

		c[1] = x(0, 0) * std::pow(x(0, 3), 3.0) / 6.0;
		c[2] = 2.0 * x(0, 0) * x(0, 3);
		c[3] = (x(0, 1) - x(0, 3) / 2.0);
		c[4] = x(0, 2) * std::pow(x(0, 1) - 2.0 * x(0, 3), 3.0) / 12.0;

		c[0] = 2.0 * x(0, 0) * x(0, 2) + x(0, 2) * (x(0, 1) - 2.0 * x(0, 3));
		for (size_t i = 0; i < 1; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 3:
		//weldead beam design 
		temp = 1.10471 * std::pow(x(0, 0), 2.0) * x(0, 1) + 0.04811 * x(0, 2) * x(0, 3) * (14 + x(0, 1));
		M = P * (L + x(0, 1) / 2.0);
		R = std::sqrt(0.25 * (x(0, 1) * x(0, 1) + (x(0, 0) + std::pow(x(0, 2), 2.0))));
		J = 2 * (std::sqrt(2) * x(0, 0) * x(0, 1) * (x(0, 1) * x(0, 1) / 12 + 0.25 * (x(0, 0) + x(0, 2) * x(0, 2))));
		P_c = (4.013 * E / (6 * L * L)) * x(0, 2) * x(0, 3) * x(0, 3) * (1 - 0.25 * x(0, 2) * std::sqrt(E / G) / L);
		t1 = P / (std::sqrt(2) * x(0, 0) * x(0, 1));
		t2 = M * R / J;
		t = std::sqrt(t1 * t1 + t1 * t2 * x(0, 1) / R + t2 * t2);
		s = 6 * P * L / (x(0, 3) * x(0, 2) * x(0, 2));
		d = 4 * P * std::pow(L, 3.0) / (E * x(0, 3) * std::pow(x(0, 2), 3.0));
		c[0] = t - t_max;
		c[1] = s - s_max;;
		c[2] = x(0, 0) - x(0, 3);
		c[3] = 0.10471 * std::pow(x(0, 0), 2.0) + 0.04811 * x(0, 2) * x(0, 3) * (14.0 + x(0, 1)) - 5.0;
		c[4] = 0.125 - x(0, 0);
		c[5] = d - d_max;
		c[6] = P - P_c;
		for (size_t i = 0; i <= 6; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 4:
		//speed reducer design 
		temp = 0.7854 * x(0, 0) * std::pow(x(0, 1), 2.0) *
			(3.3333 * std::pow(x(0, 2), 2.0) + 14.9334 * x(0, 2) - 43.0934) -
			1.508 * x(0, 0) * (std::pow(x(0, 5), 2.0) + std::pow(x(0, 6), 2.0))
			+ 7.4777 * (std::pow(x(0, 5), 3.0) + std::pow(x(0, 6), 3.0)) +
			0.7854 * (x(0, 3) * std::pow(x(0, 5), 2.0) + x(0, 4) * std::pow(x(0, 6), 2.0));
		c[0] = (27.0 / (x(0, 0) * std::pow(x(0, 1), 2.0) * x(0, 2))) - 1.0;
		c[1] = (397.5 / (x(0, 0) * std::pow(x(0, 1), 2.0) * std::pow(x(0, 2), 2.0))) - 1.0;
		c[2] = (1.93 * std::pow(x(0, 3), 4.0) / (x(0, 1) * std::pow(x(0, 5), 4.0) * x(0, 2))) - 1.0;
		c[3] = (1.93 * std::pow(x(0, 4), 3.0) / (x(0, 1) * std::pow(x(0, 6), 4.0) * x(0, 2))) - 1.0;

		c[4] = 1.0 / 110 * std::pow(x(0, 5), 3.0) * std::sqrt(std::pow(745 * x(0, 3) / x(0, 1) * x(0, 2), 2.0) + 16.9 * 1e6) - 1.0;
		c[5] = 1.0 / 85.0 * std::pow(x(0, 6), 3.0) * std::sqrt(std::pow(745 * x(0, 4) / x(0, 1) * x(0, 2), 2.0) +
			157.5 * 1e6) - 1.0;
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
	case 5:
		//cantilever beam 
		temp = 0.0;
		for (size_t i = 0; i < 5; i++)
		{
			temp += x(0, i);
		}
		temp = temp * 0.0624;
		c[0] = (61.0 / std::pow(x(0, 0), 3.0)) + (37.0 / std::pow(x(0, 1), 2.0)) +
			(19.0 / std::pow(x(0, 2), 3.0)) + (7.0 / std::pow(x(0, 3), 4.0)) +
			(1.0 / std::pow(x(0, 4), 3.0)) - 1.0;
		return (temp + PENALTY * std::pow(std::max(c[0], 0.0), 2.0));
		break;
	case 6:
		//currugated bulkhead design 
		c[0] = 5.885 * x(0, 3) * (x(0, 0) + x(0, 2));
		c[1] = std::fabs(std::pow(x(0, 2), 2.0) - std::pow(x(0, 1), 2.0));
		c[2] = x(0, 0) + std::sqrt(c[1]);
		temp = c[0] / c[2];
		//constrains
		c[3] = -1.0 * x(0, 3) * x(0, 1) * (0.4 * x(0, 0) + x(0, 2) / 6.0) + 8.94 * (x(0, 0) + std::sqrt(c[1]));
		c[4] = -1.0 * x(0, 3) * std::pow(x(0, 1), 2.0) * (0.2 * x(0, 0) + x(0, 2) / 12.0) + 2.2 * std::pow((8.94 * (x(0, 0) + std::sqrt(c[1]))), 1.33333333333);
		c[5] = -1.0 * x(0, 3) + 0.0156 * x(0, 0) + 0.15;
		c[6] = -1.0 * x(0, 3) + 0.0156 * x(0, 2) + 0.15;
		c[7] = -1.0 * x(0, 3) + 1.05;
		c[8] = -1.0 * x(0, 2) + x(0, 1);
		for (size_t i = 3; i < 9; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 7:
		//piston lever
		temp = 9.8 * x(0, 0) * x(0, 1) + 2.0 * x(0, 0);
		c[0] = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0);

		c[1] = (2500 / M_PI * x(0, 0) * x(0, 1) * 500.0) - 1.0;
		c[2] = (8.0 * 2500 * std::pow(250.0, 2.0) / (std::pow(M_PI, 3.0) * 0.85 * 106 * x(0, 0) * x(0, 1) * c[0])) - 1.0;
		c[3] = 2.0 / x(0, 0) - 1.0;
		c[4] = x(0, 0) / 14.0 - 1.0;
		c[5] = 0.2 / x(0, 1) - 1.0;
		c[6] = x(0, 1) / 0.8 - 1.0;
		for (size_t i = 1; i < 7; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 8:
		//gear trains
		c[0] = 1 / 6.931;
		c[1] = x(0, 1) * x(0, 2) / x(0, 0) * x(0, 3);
		temp = std::pow(c[0] - c[1], 2.0);
		return temp;
		break;
	case 9:
		//multiple disc clutch brake problem 
		mC = { 20.0,55.0,7.8e-6,1.0,1000.0,15.0,0.5,1.5,40.0,3.0,250.0,10.0,30.0,0.5 };
		//delr 0 ,iz 1,pau 2,pmax 3, fmax 4, tmax 5, miu 6, s 7, ms 8, mf 9, n   10, vsmax 11, lmax 12, del 13
		c[0] = (std::pow(x(0, 1), 3.0) - std::pow(x(0, 0), 3.0)); //cube term
		c[1] = (std::pow(x(0, 1), 2.0) - std::pow(x(0, 0), 2.0)); //square term 

		c[2] = (2 / 3.0) * mC[6] * x(0, 3) * x(0, 4) * c[0] / c[1];//Mh
		c[3] = x(0, 3) / M_PI * c[1];//prz
		c[4] = 2.0 * M_PI * mC[10] * c[0] / 90.0 * c[1];
		c[5] = mC[1] * M_PI * mC[10] / 30.0 * (c[2] / 1000.0 + mC[9]);
		//penalty term 
		c[6] = -mC[12] + (x(0, 4) + 1.0) * (x(0, 2) + mC[13]);
		c[7] = -mC[3] + c[3];
		c[8] = -mC[3] * mC[11] * 1000.0 + c[3] * c[4];
		c[9] = -mC[11] * 1000.0 + mC[11];
		c[10] = -mC[5] + c[5];
		c[11] = -c[2] + mC[7] * mC[8] * 1000.0;
		c[12] = -c[5];
		c[13] = -x(0, 1) + x(0, 0) + mC[0];
		temp = M_PI * c[1] * x(0, 2) * (x(0, 4) + 1.0) * mC[2];

		for (size_t i = 6; i < 14; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp+sum);
		break;
	case 10:
		//three-bar truss design 
		temp = (2.0 * std::sqrt(2.0) * x(0, 0) + x(0, 1)) * 100.0;
		c[0] = ((std::sqrt(2.0) * x(0, 0) + x(0, 1)) / (std::sqrt(2.0) * std::pow(x(0, 0), 2.0) + 2.0 * x(0, 0) * x(0, 1))) * 2.0 - 2.0;
		c[1] = ((x(0, 1)) / (std::sqrt(2.0) * std::pow(x(0, 0), 2.0) + 2.0 * x(0, 0) * x(0, 1))) * 2.0 - 2.0;
		c[2] = (1.0 / (std::sqrt(2.0) * x(0, 1) + x(0, 0))) * 2.0 - 2.0;
		for (size_t i = 0; i <= 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return temp + sum;
		break;
	case 12://problem A.1 ali sadolah.fix
		temp = std::pow(x(0, 0) - 2.0,2.0) + std::pow(x(0, 1) - 1.0, 2.0);
		c[0] = x(0, 0) - 2.0 * x(0, 1) + 1.0;//equality 
		c[1] = (std::pow(x(0, 0), 2.0) / 4.0) + std::pow(x(0, 1), 2.0) - 1.0;
		return (temp + PENALTY * std::pow(std::max(c[1], 0.0), 2.0) + PENALTY*std::pow(c[0], 2.0));
		break;
	case 13://problem a3. .fix 
		temp = std::pow(std::pow(x(0, 0), 2.0) + x(0, 1) - 11.0, 2.0) 
			+  std::pow(x(0,0)+std::pow(x(0,1),2.0)-7.0,2.0);
		c[0] = -4.84 + std::pow(x(0, 0) - 0.05, 2.0) + std::pow(x(0, 1) - 2.5, 2.0);
		c[1] = -std::pow(x(0, 0), 2.0) - std::pow(x(0, 1) - 2.5, 2.0) + 4.84;
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	case 14://problem a5 .fix
		temp = std::pow(x(0, 0) - 10.0, 3.0) + std::pow(x(0, 1) - 20.0, 3.0);
		c[0] = - std::pow(x(0, 0) - 5.0, 2.0) - std::pow(x(0, 1) - 5.0, 2.0) + 100.0;
		c[1] = std::pow(x(0, 0) - 6.0, 2.0) + std::pow(x(0, 1) - 5.0, 2.0) - 82.81;
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	
	case 15://Quantum seeded evolutionary. problem g10.fix
		temp = x(0, 0) + x(0, 1) + x(0, 2);
		c[0] = -1.0 + 0.0025 * (x(0, 3) + x(0, 5));
		c[1] = -1.0 + 0.0025 * (x(0, 4) + x(0, 6) - x(0, 3));
		c[2] = -1 + 0.01 * (x(0, 7) - x(0, 4));
		c[3] = -x(0, 0) * x(0, 5) + 833.33252 * x(0, 3) + 100 * x(0, 0) - 83333.333;
		c[4] = -x(0, 1) * x(0, 6) + 1250 * x(0, 4) + x(0, 1) * x(0, 3) - 1250 * x(0, 3);
		c[5] = -x(0, 2) * x(0, 7) + 1250000 + x(0, 2) * x(0, 4) - 2500 * x(0, 4);
		for (size_t i = 0; i < 6; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	case 16: //problem a7. sadollah stop problem  set .fix 
		temp = 5.3578547 * std::pow(x(0, 2), 2.0) + 0.8356891 * x(0, 0) * x(0, 4) +
			37.293239 * x(0, 0) - 40729.141;
		c[0] = 85.334407 + 0.0056858 * x(0, 1) * x(0, 4) + 0.0006262 * x(0, 0) * x(0, 3)
			- 0.0022053 * x(0, 2) * x(0, 4) - 92.0;
		c[1]  = -85.334407 - 0.0056858 * x(0, 1) * x(0, 4) - 0.0006262 * x(0, 0) * x(0, 3)
			-0.0022053 * x(0, 2) * x(0, 4) + 92.0;
		c[2] = 80.51249 + 0.0071317 * x(0, 1) * x(0, 4) + 0.0029955 * x(0, 0) * x(0, 1)
			+ 0.0021813 * std::pow(x(0, 2), 2.0) - 110.0;
		c[3] = -80.51249 - 0.0071317 * x(0, 1) * x(0, 4) - 0.0029955 * x(0, 0) * x(0, 1)
			- 0.0021813 * std::pow(x(0, 2), 2.0) + 90.0;
		c[4] = 9.300961 + 0.0047026 * x(0, 2) * x(0, 4) + 0.0012547*x(0, 0) * x(0, 2)
		+0.0019085*x(0, 2) * x(0, 3) - 25.0;
		c[5] = -9.300961 - 0.0047026 * x(0, 2) * x(0, 4) - 0.0012547 * x(0, 0) * x(0, 2)
		-0.0019085 * x(0, 2) * x(0, 3) + 20.0;
		for (size_t i = 0; i < 6; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	case 17://qiang zhao, two stage multi swarm.fix problem 5 fix 
		temp = std::pow(x(0, 0) - 10.0, 2.0) + 5.0 * std::pow(x(0, 1) - 12.0, 2.0) + std::pow(x(0, 2), 4.0) +
			3.0 * std::pow(x(0, 3) - 11.0, 2.0) + 10.0 * std::pow(x(0, 4), 6.0) + 7.0 * std::pow(x(0, 5), 2.0) +
			std::pow(x(0, 6), 4.0) - 4.0 * x(0, 5) * x(0, 6) - 10.0 * x(0, 5) - 8.0 * x(0, 6);
		c[1] = -127.0 + 2.0 * std::pow(x(0, 0), 2.0) + 3.0 * std::pow(x(0, 1), 4.0) + x(0, 2) + 4.0 * std::pow(x(0, 3), 2.0) +
			5.0 * x(0, 4);
		c[2] = -282.0 + 7.0 * x(0, 0) + 3.0 * x(0, 1) + 10.0 * std::pow(x(0, 2), 2.0) + x(0, 3) - x(0, 4);
		c[3] = -196.0 + 23.0 * x(0, 0) + std::pow(x(0, 1), 2.0) + 6.0 * std::pow(x(0, 5), 2.0) - 8.0 * x(0, 6);
		c[4] = 4.0 * std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) - 3.0 * x(0, 0) * x(0, 1) +
			2.0 * std::pow(x(0, 2), 2.0) + 5.0 * x(0, 5) - 11.0 * x(0, 6);
		for (size_t i = 1; i < 5; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (temp + sum);
		break;
	case 18://qiang zhao, two stage multi swarm.fix
		//b03
		temp = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) + x(0, 0) * x(0, 1) - 14.0 * x(0, 0) - 16.0 * x(0, 1) +
			std::pow(x(0, 2) - 10.0, 2.0) + 4.0 * std::pow(x(0, 3) - 5.0, 2.0) + std::pow(x(0, 4) - 3.0, 2.0) +
			2.0 * std::pow(x(0, 5) - 1.0, 2.0) + 5.0 * std::pow(x(0, 6), 2.0) + 7.0 * std::pow(x(0, 7) - 11.0, 2.0) +
			2.0 * std::pow(x(0, 8) - 10.0, 2.0) + std::pow(x(0, 9) - 7.0, 2.0) + 45.0;
		c[0] = -105.0 + 4.0 * x(0, 0) + 5.0 * x(0, 1) - 3.0 * x(0, 6) + 9.0 * x(0, 7);
		c[1] = 10 * x(0, 0) - 8.0 * x(0, 1) - 17.0 * x(0, 6) + 2.0 * x(0, 7);
		c[2] = -8.0 * x(0, 0) + 2.0 * x(0, 1) + 5.0 * x(0, 8) - 2.0 * x(0, 9) - 12.0;
		c[3] = 3.0 * std::pow(x(0, 0) - 2.0, 2.0) + 4.0 * std::pow(x(0, 1) - 3.0, 2.0) + 2.0 * std::pow(x(0, 2), 2.0) - 7.0 * x(0, 3)-120.0;
		c[4] = 5.0 * std::pow(x(0, 0), 2.0) + 8.0 * x(0, 1) + std::pow(x(0, 2) - 6.0, 2.0) - 2.0 * x(0, 3) - 40.0;
		c[5] = std::pow(x(0, 0), 2.0) + 2.0 * std::pow(x(0, 1) - 2.0, 2.0) - 2.0 * x(0, 0) * x(0, 1) + 14.0 * x(0, 4) - 6.0 * x(0, 5);
		c[6] = 0.5 * std::pow(x(0, 0) - 8.0, 2.0) + 2.0 * std::pow(x(0, 1) - 4.0, 2.0) + 3.0 * std::pow(x(0, 4), 2.0) - x(0, 5) - 30.0;
		c[7] = -3.0 * x(0, 0) + 6.0 * x(0, 1) + 12.0 * std::pow(x(0, 8) - 8.0, 2.0) - 7.0 * x(0, 9);
		for (size_t i = 0; i < 8; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	case 19://b01 qiang zhao, two stage multi swarm ..fix
		temp = 3.0 * x(0, 0) + 0.000001 * std::pow(x(0, 0), 2.0) + 2.0 * x(0, 1) 
			+ (0.000002 / 3.0) * std::pow(x(0, 1), 3.0);
		c[0] = -x(0, 3) + x(0, 2) - 0.55; //inequality 
		c[1] = -x(0, 2) + x(0, 3) - 0.55;
		c[2] = 1000.0 * std::sin(-x(0, 2) - 0.25) + 1000.0 * std::sin(-x(0, 3) - 0.25) + 894.8 - x(0, 0);
		c[3] = 1000.0 * std::sin(x(0, 2) - 0.25) + 1000.0 * std::sin(x(0, 2)-x(0,3) - 0.25) + 894.8 - x(0, 1);
		c[4] = 1000.0 * std::sin(x(0, 3) - 0.25) + 1000.0 * std::sin(x(0, 3) - x(0, 2) - 0.25) + 1294.8;
		c[9] = 0.0;
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		for (size_t i = 2; i < 5; i++)
		{
			c[9] += PENALTY * std::pow(c[i],2.0);
		}
		return	(temp + sum+c[9]);
		break;
	case 20://A novel filled function method for solving.qiang li. fix. problem 4.1  
		temp = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) - std::cos(17.0 * x(0, 0)) - std::cos(17.0 * x(0, 1)) + 3.0;
		c[0] = std::pow(x(0, 0) - 2.0, 2.0) + std::pow(x(0, 1), 2.0) - std::pow(1.6, 2.0);
		c[1] = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1) - 3.0, 2.0) - std::pow(2.7, 2.0);
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	case 21://A novel filled function method for solving.qiang li.  problem 4.2   fix.
		temp = -x(0, 0) - x(0, 1);
		c[0] = x(0, 1) - 2.0 * std::pow(x(0, 0), 4.0)
			+ 8.0 * std::pow(x(0, 0), 3.0) - 8.0 * std::pow(x(0, 0), 2.0) - 2.0;
		c[1] = x(0, 1) - 4.0 * std::pow(x(0, 0), 4.0)
			+ 32 * std::pow(x(0, 0), 3.0) - 88.0 * std::pow(x(0, 0), 2.0) + 96.0 * x(0, 0) - 36.0;
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	case 22: // Value - Estimation Function - x.l sun.problem 2.fix
		temp = std::pow(x(0, 0) - 2.0, 2.0) + std::pow(x(0, 1) - 1.0, 2.0) +1.0;
		c[0] = x(0, 0) + x(0, 1) - 2.0;
		c[1] = std::pow(x(0, 0), 2.0) - x(0, 1);
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	case 23://A novel modified BSA, hailong wang, problem 3.fix
		temp = std::pow(std::pow(x(0, 0), 2.0) + x(0, 1) - 11.0, 2.0) + std::pow(std::pow(x(0, 1), 2.0) + x(0, 0) - 7.0, 2.0);
		c[0] = std::pow(x(0, 0) - 0.05, 2.0) + std::pow(x(0, 1) - 2.5, 2.0) - 4.84;
		c[1] = -std::pow(x(0, 0), 2.0) - std::pow(x(0, 1) - 2.5, 2.0) + 4.84;
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break; 
	case 24://A novel modified BSA, hailong wang, problem 4.fix
		temp = -std::pow(std::sin(2.0 * M_PI * x(0, 0)), 3.0) * (std::sin(2.0 * M_PI * x(0, 1)))
			/ std::pow(x(0, 0), 3.0) * (x(0, 0) + x(0, 1));
		c[0] = std::pow(x(0, 0), 2.0) - x(0, 1) + 1.0;
		c[1] = 1.0 - x(0, 0) + std::pow(x(0, 1) - 4.0, 2.0);
		for (size_t i = 0; i < 2; i++)
		{
			sum += PENALTY * std::pow(std::max(c[i], 0.0), 2.0);
		}
		return (sum + temp);
		break;
	
	/*case 0:
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
	
	//multiobjective optimization 
	case 35:
		c[0] = std::pow(x(0, 0),2.0);
		c[1] = std::pow(x(0, 0) - 2.0, 2.0);
		for (size_t i = 0; i < 2; i++)
		{
			sum += 0.5 * c[i];
		}
		return sum;
		break;
		*/
	case 25: //schewefel function version 2 
		temp = 0;
		for (size_t i = 0; i < nvars; i++)
		{
			c[0] = x(0, i);
			temp = temp + c[0]*std::sin(std::sqrt(std::fabs(c[0])));
		}
		return (418.9829 * nvars - temp);
		break;
	case 26://rastrign function  
		for (size_t i = 0; i < nvars; i++)
		{
			c[0] = x(0, i);
			temp += std::pow(c[0], 2.0) - 10 * std::cos(2.0 * M_PI * c[0]) + 10.0;
		}
		return temp;
		break;
	case 27://ackley 
		c[0] = 0, c[1] = 0,
			c[2]=0,c[3]=0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp = x(0, i);
			c[0] += std::pow(temp, 2.0);
			c[1] += std::cos(2.0 * M_PI * temp);
		}
		c[2] = -20 * std::exp(-0.2 * std::sqrt(c[0] / nvars));
		c[3] = -std::exp(c[1] / nvars);
		return (c[2] + c[3] + 20 + std::exp(1.0));
		break;
	case 28:
		c[0] = 0, c[1] = 1.0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp = x(0, i);
			c[0] += std::pow(temp, 2.0) / 4000;
			c[1] *= std::cos(temp / std::sqrt(i + 1));
		}
		return (c[0] - c[1] + 1.0);
		break;
		//bencmark 2 dimension constraint
	case 100://rosenbrock 
		temp = std::pow(1 - x(0, 0), 2.0) + 100 * std::pow(x(0, 1) - x(0, 0) * x(0, 0), 2.0);
		c[0] = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) - 2.0;
		return (temp + beta * std::pow(std::max(0.0, c[0]),2.0));
		break;
	case 30://rosenbrock 
		temp = std::pow(1 - x(0, 0), 2.0) + 100 * std::pow(x(0, 1) - x(0, 0) * x(0, 0), 2.0);
		c[0] = std::pow(x(0, 0) - 1.0, 3.0) - x(0, 1) + 1.0;
		c[1] = x(0, 0) + x(0, 1) - 2.0;
		c[2] = 0.0;
		for (size_t i = 0; i < 2; i++)
		{
			c[2] += beta * std::pow(std::max(0.0, c[i]), 2.0);
		}
		return (temp + c[2]);
		break;
	case 31://simionescu
		temp = 0.1 * x(0, 0) * x(0, 1);
		c[0] = std::pow( 1.0 + 0.2 * std::cos(8.0 * std::atan(x(0, 0) / x(0, 1))),2.0);
		c[1] = std::pow(x(0, 0), 2.0) + std::pow(x(0, 1), 2.0) - c[0];
		return (temp + beta * std::pow(std::max(0.0, c[1]), 2.0));
		break;
	case 29://mishra's Bird
		c[0] = std::pow(1.0 - std::cos(x(0, 0)), 2.0);
		c[1] = std::pow(1.0 - std::sin(x(0, 1)), 2.0);
		temp = std::sin(x(0, 1)) * std::exp(c[0]) + 
			std::cos(x(0, 0)) * std::exp(c[1]) + std::pow(x(0, 0) - x(0, 1), 2.0);
		c[2] = std::pow(x(0, 0) + 5.0, 2.0) + std::pow(x(0, 1) + 5.0, 2.0) - 25.0;
		return (temp + beta * std::pow(std::max(0.0, c[2]), 2.0));
		break;
	default:
		break;
	}
	return 0;
}
void NNA::initialization(Eigen::Matrix<double, 1, NPOP>& ww, Eigen::MatrixXd& w,
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
			switch (KODE)
			{
				//tension 
			case 0:
				if ((j == 0))
				{
					x_pattern(i, j) = tn1(eng);
				}
				else if (j == 1) {
					x_pattern(i, j) = tn2(eng);
				}
				else {
					x_pattern(i, j) = tn3(eng);
				}
				break;
			case 1:
				//pressure vessel
				if ((j == 0) || (j == 1))
				{
					x_pattern(i, j) = pr1(eng);
				}
				else {
					x_pattern(i, j) = pr3(eng);
					}
				break;
			case 2:
				// i-beam design 
				if ((j == 2) || (j == 3))
				{
					x_pattern(i, j) = ib3(eng);
				}
				else if (j==1)
				{
					x_pattern(i, j) = ib2(eng);
				}
				else {
					x_pattern(i, j) = ib1(eng);
				}
				break;
			case 3:
				//weldead beam
				if ((j == 0) || (j == 3))
				{
					x_pattern(i, j) = wb1(eng);
				}
				else {
					x_pattern(i, j) = wb2(eng);

				}
				break;
			case 4:
				//speed reducer 
				if (j==0)
				{
					x_pattern(i, j) = d1(eng);
				}
				else if (j == 1)
				{
					x_pattern(i, j) = d2(eng);
				}
				else if (j == 2)
				{
					x_pattern(i, j) = d3(eng);
				}
				else if (j == 3)
				{
					x_pattern(i, j) = d4(eng);
				}
				else if (j == 4)
				{
					x_pattern(i, j) = d5(eng);
				}
				else if (j == 5)
				{
					x_pattern(i, j) = d6(eng);
				}
				else {
					x_pattern(i, j) = d7(eng);
				}
				break;
			case 5:
				//cantvilear beam design 
				x_pattern(i, j) = cb(eng); 
				break;
			case 6:
				//corrugated bulkhead design 
				if ( (j == 0) || (j == 1) || (j == 2) )
				{
					x_pattern(i, j) = cur1(eng);
				}
				else
				{
					x_pattern(i, j) = cur2(eng);
				}
				break;
			case 7:
				//piston lever
				x_pattern(i, j) = tC(eng); 
				break;
			case 8:
				//gear train 
				x_pattern(i, j) = ge(eng);
				x_pattern(i, j) = std::round(x_pattern(i, j)); 
				break;
			case 9:
				//multiple disk 
				if (j == 0)
				{
					x_pattern(i, j) = md1(eng);
					x_pattern(i, j) = std::round(x_pattern(i, j));
				}
				else if (j == 1)
				{
					x_pattern(i, j) = md2(eng);
					x_pattern(i, j) = std::round(x_pattern(i, j));

				}
				else if (j == 2)
				{
					x_pattern(i, j) = bT[muldiskT(eng)];

				}
				else if (j == 3)
				{
					x_pattern(i, j) = bF[muldiskF(eng)];

				}
				else
				{
					x_pattern(i, j) = md5(eng);
					x_pattern(i, j) = std::round(x_pattern(i, j));

				}
				break;
			case 14:
				if (j==0)
				{
					x_pattern(i, j) = CS1(eng);
				}
				else
				{
					x_pattern(i, j) = CS2(eng);
				}
				break;
			case 15:
				if (j == 0)
				{
					x_pattern(i, j) = DS1(eng);
				}
				else if (j == 1|| j==2)
				{
					x_pattern(i, j) = DS2(eng);
				}
				else
				{
					x_pattern(i, j) = DS2(eng);
				}
				break;
			case 16:
				if (j == 0)
				{
					x_pattern(i, j) = ES1(eng);
				}
				else if (j == 1)
				{
					x_pattern(i, j) = ES2(eng);
				}
				else
				{
					x_pattern(i, j) = ES3(eng);
				}
				break;
			case 19:
				if (j == 0)
				{
					x_pattern(i, j) = FS1(eng);
				}
				else if (j == 1)
				{
					x_pattern(i, j) = FS2(eng);
				}
				else
				{
					x_pattern(i, j) = FS3(eng);
				}
				break;
			case 21:
				if (j == 0)
				{
					x_pattern(i, j) = GS1(eng);
				}
				else
				{
					x_pattern(i, j) = GS2(eng);
				}
				break;
			case 29://case 29-32 constrained 2-dimension 
				if (j == 0)
				{
					x_pattern(i, j) = CS29_32(eng);
				}
				else
				{
					x_pattern(i, j) = DS29_32(eng);
				}
				break;
			default:	
				x_pattern(i, j) = DF(eng);
				break;
			}
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
void NNA::generateWeight(Eigen::MatrixXd& w, Eigen::Matrix<double, 1, NPOP - 1>& t) {
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
	Eigen::MatrixXd& w) {
	x_target = x_pattern.row(int(value_index(0, 1)));
	target = value_index(0, 0);
	w_target = w.col(int(value_index(0, 1)));
}
void NNA::Run(Eigen::Matrix<double, NPOP, nvars>& x_new, Eigen::Matrix<double, NPOP, nvars>& x_pattern, Eigen::Matrix<double, NPOP, 1>& w_target,
	Eigen::MatrixXd& w, Eigen::Matrix<double, NPOP, 1>& cost, Eigen::Matrix<double, 1, nvars>& x_target) {
	std::ofstream writeF;
	int Maximum = 151;
	int begin = 1;
	int auxilaryVar = 0;
	writeF.open("case29_32.csv", std::ios::app);
	while (begin != MAX)
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
					auxilaryVar = newDistr(eng);
					//modify here
					switch (KODE)
					{
					case 0:
						//tension
						if (2 == auxilaryVar)
						{
							x_pattern(i, auxilaryVar) = tn3(eng);
						}
						else if (1== auxilaryVar)
						{
							x_pattern(i, auxilaryVar) = tn2(eng);
						}
						else {
							x_pattern(i, auxilaryVar) = tn1(eng);
						}
						break;
					case 1:
						//pressure Vessel
						if ((auxilaryVar==1)|| (auxilaryVar==2))
						{
							//x_pattern(i, auxilaryVar) = L3 + (U3 - L3) * distr(eng);
							x_pattern(i, auxilaryVar) = pr1(eng);
						}
						else {
							x_pattern(i, auxilaryVar) = pr2(eng);
						}
						break;
					case 2:
						//i-Beam 
						if ((auxilaryVar == 2) || (auxilaryVar == 3))
						{
							x_pattern(i, auxilaryVar) = ib3(eng);
						}
						else if (auxilaryVar==1)
						{
							x_pattern(i, auxilaryVar) = ib2(eng);
						}
						else {
							x_pattern(i, auxilaryVar) = ib1(eng);
						}
						break;
					case 3:
						//weldead beam 
						if ((auxilaryVar == 0) || (auxilaryVar == 3))
						{
							x_pattern(i, auxilaryVar) = wb1(eng);
						}
						else {
							x_pattern(i, auxilaryVar) = wb2(eng);
						} 
						break;
					case 4:
						//speed reducer 
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = d1(eng);
						}
						else if (auxilaryVar == 1)
						{
							x_pattern(i, auxilaryVar) = d2(eng);
						}
						else if (auxilaryVar == 2)
						{
							x_pattern(i, auxilaryVar) = d3(eng);
						}
						else if (auxilaryVar == 3)
						{
							x_pattern(i, auxilaryVar) = d4(eng);
						}
						else if (auxilaryVar == 4)
						{
							x_pattern(i, auxilaryVar) = d5(eng);
						}
						else if (auxilaryVar == 5)
						{
							x_pattern(i, auxilaryVar) = d6(eng);
						}
						else {
							x_pattern(i, auxilaryVar) = d7(eng);
						}
						break;
					case 5:
						//cantevelar beam 
						x_pattern(i, auxilaryVar) = cb(eng); 
						break;
					case 6:
						//corrugated buckhead 
						if ((auxilaryVar == 0) || (auxilaryVar == 1) || (auxilaryVar == 2))
						{
							x_pattern(i, auxilaryVar) = cur1(eng);
						}
						else {
							x_pattern(i, auxilaryVar) = cur2(eng);
						} 
						break;
					case 7:
						//piston lever
						x_pattern(i, auxilaryVar) = tC(eng); 
						break;
					case 8:
						//gear train 
						x_pattern(i, auxilaryVar) = ge(eng);
						x_pattern(i, auxilaryVar) = std::round(x_pattern(i, auxilaryVar)); 
						break;
					case 9:
						//multiple disk
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = md1(eng);
							x_pattern(i, auxilaryVar) = std::round(x_pattern(i, auxilaryVar));

						}
						else if (auxilaryVar == 1)
						{
							x_pattern(i, auxilaryVar) = md2(eng);
							x_pattern(i, auxilaryVar) = std::round(x_pattern(i, auxilaryVar));

						}
						else if (auxilaryVar == 2)
						{
							x_pattern(i, auxilaryVar) = bT[muldiskT(eng)];

						}
						else if (auxilaryVar == 3)
						{
							x_pattern(i, auxilaryVar) = bF[muldiskF(eng)];
							x_pattern(i, auxilaryVar) = std::round(x_pattern(i, auxilaryVar));

						}
						else
						{
							x_pattern(i, auxilaryVar) = md5(eng);
							x_pattern(i, auxilaryVar) = std::round(x_pattern(i, auxilaryVar));

						}
						break;
					case 14:
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = CS1(eng);
						}
						else  
						{
							x_pattern(i, auxilaryVar) = CS2(eng);
						}
						break;
					case 15:
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = DS1(eng);
						}
						else if (auxilaryVar == 1 || auxilaryVar==2)
						{
							x_pattern(i, auxilaryVar) = DS2(eng);
						}
						else
						{
							x_pattern(i, auxilaryVar) = DS3(eng);
						}
						break;
					case 16:
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = ES1(eng);
						}
						else if (auxilaryVar == 1)
						{
							x_pattern(i, auxilaryVar) = ES2(eng);
						}
						else
						{
							x_pattern(i, auxilaryVar) = ES3(eng);
						}
						break;
					case 19:
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = FS1(eng);
						}
						else if (auxilaryVar == 1)
						{
							x_pattern(i, auxilaryVar) = FS2(eng);
						}
						else
						{
							x_pattern(i, auxilaryVar) = FS3(eng);
						}
						break;
					case 21:
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = GS1(eng);
						}
						else 
						{
							x_pattern(i, auxilaryVar) = GS2(eng);
						}
						break;
					case 29://case 29-32-constrained oprimization
						if (auxilaryVar == 0)
						{
							x_pattern(i, auxilaryVar) = CS29_32(eng);
						}
						else
						{
							x_pattern(i, auxilaryVar) = DS29_32(eng);
						}
						break;
					default:
						x_pattern(i, auxilaryVar) = DF(eng);
						break;
					}
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
					//x_pattern(i, p) = std::round(x_pattern(i, p)); //mix integer-gear train,multiple disk
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
				switch (KODE)
				{
				case 0:
					//tension ...
					if (y_max == 0)
					{
						if ( (x_pattern(x_max, y_max) < tension[0] ) || (x_pattern(x_max, y_max) > tension[1]) )
						{
							x_pattern(x_max, y_max) = tn1(eng);
						}
					}
					else if (y_max == 1) {
						if ((x_pattern(x_max, y_max) < tension[2]) || (x_pattern(x_max, y_max) > tension[3]))
						{

							x_pattern(x_max, y_max) = tn2(eng);
						}
					}
					else {
						if ((x_pattern(x_max, y_max) < tension[4]) || (x_pattern(x_max, y_max) > tension[5]))
						{
							x_pattern(x_max, y_max) = tn3(eng);
						}
					}
					break;
				case 1:
					//pressure Vessel design 
					if ((y_max == 0) || (y_max==1))
					{
						if ((x_pattern(x_max, y_max) < pressureVes[0]) || (x_pattern(x_max, y_max) > pressureVes[1]))
						{
							x_pattern(x_max, y_max) = pr1(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < pressureVes[4]) || (x_pattern(x_max, y_max) > pressureVes[5]))
						{
							x_pattern(x_max, y_max) = pr3(eng);
						}
					}
					break;
				case 2:
					//ibeam design 
					if ((y_max == 2) || (y_max == 3))
					{
						if ((x_pattern(x_max, y_max) < ibeam[4]) || (x_pattern(x_max, y_max) > ibeam[5]))
						{
							x_pattern(x_max, y_max) = ib3(eng);
						}
					}
					else if (y_max == 1)
					{
						if ((x_pattern(x_max, y_max) < ibeam[2]) || (x_pattern(x_max, y_max) > ibeam[3]))
						{
							x_pattern(x_max, y_max) = ib2(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < ibeam[0]) || (x_pattern(x_max, y_max) > ibeam[1]))
						{
							x_pattern(x_max, y_max) = ib1(eng);
						}
					}
					break;
				case 3:
					//weldead beam 
					if ((y_max == 0) || (y_max == 3))
					{
						if ((x_pattern(x_max, y_max) < weldeadBeam[0]) || (x_pattern(x_max, y_max) > weldeadBeam[1]))
						{
							x_pattern(x_max, y_max) = wb1(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < weldeadBeam[2]) || (x_pattern(x_max, y_max) > weldeadBeam[3]))
						{
							x_pattern(x_max, y_max) = wb2(eng);
						}
					} 
					break;
				case 4:
					//speed reducer
					if (y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < LoSpeedRed[0]) ||
							(x_pattern(x_max, y_max) > UpSpeedRed[0]))
						{
							x_pattern(x_max, y_max) = d1(eng);
						}
					}
					else if (y_max == 1)
					{
						if ((x_pattern(x_max, y_max) < LoSpeedRed[1]) ||
							(x_pattern(x_max, y_max) > UpSpeedRed[1]))
						{
							x_pattern(x_max, y_max) = d2(eng);
						}
					}
					else if (y_max == 2)
					{
						if ((x_pattern(x_max, y_max) < LoSpeedRed[2]) ||
							(x_pattern(x_max, y_max) > UpSpeedRed[2]))
						{
							x_pattern(x_max, y_max) = d3(eng);
						}
					}
					else if (y_max == 3)
					{
						if ((x_pattern(x_max, y_max) < LoSpeedRed[3]) ||
							(x_pattern(x_max, y_max) > UpSpeedRed[3]))
						{
							x_pattern(x_max, y_max) = d4(eng);
						}
					}
					else if (y_max == 4)
					{
						if ((x_pattern(x_max, y_max) < LoSpeedRed[4]) ||
							(x_pattern(x_max, y_max) > UpSpeedRed[4]))
						{
							x_pattern(x_max, y_max) = d5(eng);
						}
					}
					else if (y_max == 5)
					{
						if ((x_pattern(x_max, y_max) < LoSpeedRed[5]) ||
							(x_pattern(x_max, y_max) > UpSpeedRed[5]))
						{
							x_pattern(x_max, y_max) = d6(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < LoSpeedRed[6]) ||
							(x_pattern(x_max, y_max) > UpSpeedRed[6]))
						{
							x_pattern(x_max, y_max) = d7(eng);
						}
					}
					break;
				case 5:
					//cantevelar beam 
					if ((x_pattern(x_max, y_max) < cantBeam[0]) ||
						(x_pattern(x_max, y_max) > cantBeam[1]))
					{
						x_pattern(x_max, y_max) = cb(eng);
					} 
					break;
				case 6:
					//corrugated bulkhead 
					 if ((y_max == 0) || (y_max == 1) || (y_max == 2))
						{
							if ((x_pattern(x_max, y_max) < currugatedBulkhead[0]) ||
								(x_pattern(x_max, y_max) > currugatedBulkhead[1]))
							{
								x_pattern(x_max, y_max) = cur1(eng);
							}
						}
					 else
					 {
						 if ((x_pattern(x_max, y_max) < currugatedBulkhead[2]) ||
							 (x_pattern(x_max, y_max) > currugatedBulkhead[3]))
						 {
							 x_pattern(x_max, y_max) = cur2(eng);
						 }
					 } 
					break;
				case 7:
					//piston leaver 
					if ((x_pattern(x_max, y_max) < tabColumn[0]) ||
						(x_pattern(x_max, y_max) > tabColumn[1]))
					{
						x_pattern(x_max, y_max) = tC(eng);
					}
					break;
				case 8:

					//gear train
					if ((x_pattern(x_max, y_max) < gear[0]) ||
						(x_pattern(x_max, y_max) > gear[1]))
					{
						x_pattern(x_max, y_max) = ge(eng);
						x_pattern(x_max, y_max) = std::round(x_pattern(x_max, y_max));

					}
					break;
				case 9:
					if (y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < LoMulD[0]) ||
							(x_pattern(x_max, y_max) > UpMulD[0]))
						{
							x_pattern(x_max, y_max) = md1(eng);
							x_pattern(x_max, y_max) = std::round(x_pattern(x_max, y_max));

						}
					}
					else if (y_max == 1)
					{
						if ((x_pattern(x_max, y_max) < LoMulD[1]) ||
							(x_pattern(x_max, y_max) > UpMulD[1]))
						{
							x_pattern(x_max, y_max) = md2(eng);
							x_pattern(x_max, y_max) = std::round(x_pattern(x_max, y_max));

						}
					}
					else if (y_max == 2)
					{
						if ((x_pattern(x_max, y_max) < LoMulD[2]) ||
							(x_pattern(x_max, y_max) > UpMulD[2]))
						{
							x_pattern(x_max, y_max) = bT[muldiskT(eng)];

						}
					}
					else if (y_max == 3)
					{
						if ((x_pattern(x_max, y_max) < LoMulD[3]) ||
							(x_pattern(x_max, y_max) > UpMulD[3]))
						{
							x_pattern(x_max, y_max) = bF[muldiskF(eng)];

						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < LoMulD[4]) ||
							(x_pattern(x_max, y_max) > UpMulD[4]))
						{
							x_pattern(x_max, y_max) = md5(eng);
							x_pattern(x_max, y_max) = std::round(x_pattern(x_max, y_max));

						}
					}
					break;
				case 14:
					if(y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < cs1421[0]) ||
							(x_pattern(x_max, y_max) > cs1421[1]))
						{
							x_pattern(x_max, y_max) = CS1(eng);
						}
					}
					else 
					{
						if ((x_pattern(x_max, y_max) < cs1421[2]) ||
							(x_pattern(x_max, y_max) > cs1421[3]))
						{
							x_pattern(x_max, y_max) = CS2(eng);

						}
					}
					break;
				case 15:
					if (y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < cs15[0]) ||
							(x_pattern(x_max, y_max) > cs15[1]) )
						{
							x_pattern(x_max, y_max) = DS1(eng);
						}
					}
					else if (y_max == 1 ||y_max==2)
					{
						if ((x_pattern(x_max, y_max) < cs15[2]) ||
							(x_pattern(x_max, y_max) > cs15[3]))
						{
							x_pattern(x_max, y_max) = DS2(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < cs15[4]) ||
							(x_pattern(x_max, y_max) > cs15[5]))
						{
							x_pattern(x_max, y_max) = DS3(eng);

						}
					}
					break;
				case 16:
					if (y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < cs16[0]) ||
							(x_pattern(x_max, y_max) > cs16[1]))
						{
							x_pattern(x_max, y_max) = ES1(eng);
						}
					}
					else if (y_max == 1)
					{
						if ((x_pattern(x_max, y_max) < cs16[2]) ||
							(x_pattern(x_max, y_max) > cs16[3]))
						{
							x_pattern(x_max, y_max) = ES2(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < cs16[4]) ||
							(x_pattern(x_max, y_max) > cs16[5]))
						{
							x_pattern(x_max, y_max) = ES3(eng);

						}
					}
					break;
				case 19:
					if (y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < cs19[0]) ||
							(x_pattern(x_max, y_max) > cs19[1]))
						{
							x_pattern(x_max, y_max) = FS1(eng);
						}
					}
					else if (y_max == 1)
					{
						if ((x_pattern(x_max, y_max) < cs19[2]) ||
							(x_pattern(x_max, y_max) > cs19[3]))
						{
							x_pattern(x_max, y_max) = FS2(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < cs19[4]) ||
							(x_pattern(x_max, y_max) > cs19[5]))
						{
							x_pattern(x_max, y_max) = FS3(eng);

						}
					}
					break;
				case 21:
					if (y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < cs21[0]) ||
							(x_pattern(x_max, y_max) > cs21[1]))
						{
							x_pattern(x_max, y_max) = GS1(eng);
						}
					}
					else 
					{
						if ((x_pattern(x_max, y_max) < cs21[2]) ||
							(x_pattern(x_max, y_max) > cs21[3]))
						{
							x_pattern(x_max, y_max) = GS2(eng);
						}
					}
					break;
				case 29://case 29-32 constrained optimization 
					if (y_max == 0)
					{
						if ((x_pattern(x_max, y_max) < cs29_32[0]) ||
							(x_pattern(x_max, y_max) > cs29_32[1]))
						{
							x_pattern(x_max, y_max) = CS29_32(eng);
						}
					}
					else
					{
						if ((x_pattern(x_max, y_max) < cs29_32[2]) ||
							(x_pattern(x_max, y_max) > cs29_32[3]))
						{
							x_pattern(x_max, y_max) = CS29_32(eng);
						}
					}
					break;
				default:
					if ((x_pattern(x_max, y_max) < def[0]) ||
						(x_pattern(x_max, y_max) > def[1]))
					{
						x_pattern(x_max, y_max) = DF(eng);
					}
					break;
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
		writeF << std::endl;
		//std::cout << target << std::endl;

		//std::cout << x_target << "\t" << target << std::endl;	
		++begin;
	}
	writeF <<std::endl;
	writeF.close(); 
	std::cout << x_target << std::endl;

	//std::cout << target << std::endl;
	//std::cout << "xtarge \t: " <<x_target << std::endl;
}
double fx(Eigen::Matrix<double, 1, 2>x) {
	return (std::cos(x(0, 0)) * std::exp(x(0, 1)));
	//return (std::pow(x(0, 0), 2.0) * x(0, 1));
}
