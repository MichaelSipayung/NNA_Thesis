#include <iostream>
#include <Eigen/dense>
#include <random> 
#include <algorithm>
#include "AnnImp.cpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#define NPOP 100
#define MAX  150
#define LB 0 
#define UP 10
#define nvars 2
double f(Eigen::Matrix<double, 1, nvars>);

int main()
{
	int begin = 0;
	Eigen::Matrix<double, NPOP, 1> x_lb;
	Eigen::Matrix<double, NPOP, 1> x_ub;
	Eigen::Matrix<double, NPOP, 1> cost;
	Eigen::Matrix<double, NPOP, nvars> x_pattern;
	Eigen::Matrix<double, 1, NPOP> ww; 
	Eigen::Matrix<double, NPOP, NPOP> w;
	Eigen::Matrix<double, 1, NPOP-1> t;
	Eigen::Matrix <double, 1, 2> value_index; 
	Eigen::Matrix <double, 1, 2> value_index1;
	Eigen::Matrix <double, 1, 2> value_index2;

	Eigen::Matrix<double, 1, nvars> x_target;
	Eigen::Matrix<double, NPOP, 1> w_target;
	Eigen::Matrix<double, NPOP, nvars> x_new;
	Eigen::Matrix<double, MAX, nvars> f_min;
	Eigen::Matrix<double, 1, nvars> xx;
	int  rotate_position[nvars];
	for (size_t i = 1; i < nvars+1; i++)
	{
		rotate_position[i - 1] = i;
	}

	f_min = f_min.setZero();
	int n_rotate = 0,n_wrotate=0;
	double target = 0.0; 
	double beta = 1.0; 
	constexpr int FLOAT_MIN = 0;
	constexpr int FLOAT_MAX = 1;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);
	std::uniform_int_distribution<int> newDistr(0, nvars-1), newDistrWeight(0, NPOP - 1);
	
	//step 1: initialization 
	ww = ww.setOnes()*0.5;
	w = 0.5*w.setIdentity(); 
	int var = 0;
	x_lb		= x_lb.setConstant(LB);
	x_ub		= x_ub.setConstant(UP);
	x_pattern	= x_pattern.setZero();
	cost		= cost.setZero();
	//step 1: stop, return default value. stop 
	
	//step 2 and 3 : randomly generate initial population between lower and upper bound
	// calculate the objective function
	for (size_t i = 0; i < NPOP; i++)
	{
		for (size_t j = 0; j < nvars; j++)
		{
			x_pattern(i, j) = LB + (UP - LB) * distr(eng);
		}
		cost(i,0) = f(x_pattern.row(i));
	}
	for (size_t i = 0; i < cost.rows(); i++)
	{
		if (cost(i,0)==cost.minCoeff()) {
			value_index(0, 0) = cost(i,0);
			value_index(0, 1) = i;
			break;
		}
	}
	//return the best cost and value initial population(index the best (x[1],..... x[n] appear). stop 

	//step 4 : randomly generate the weight matrix
	for (var=0; var < w.cols(); var++)
	{
		double total = 0.0;

		for (int j = 0; j < t.cols(); j++)
		{
			t(0, j) = distr(eng) * 0.5;
			total += t(0, j);

		}
		for (int k = 0; k < t.cols(); k++)
		{
			t(0, k) = (t(0,k)/total)*0.5;
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
	//step 4: return weight matrix, which  is satisfy the constrains. stop 

	//step 5: set target solution..
	x_target = x_pattern.row(int(value_index(0, 1)));
	target = value_index(0, 0);
	w_target = w.col(int(value_index(0, 1)));
	//return target_vector, target value, and weight matrix;. stop 

	while (begin!=MAX)
	{
		//step 6: generate new pattern and update solution... 
		x_new		= w * x_pattern;
		x_pattern	= x_new + x_pattern;
		//return new_solution.stop 
		
		//step 7: update the weight matrix..
		for (size_t i = 0; i < NPOP; i++)
		{
			for (size_t j = 0; j < NPOP; j++)
			{
				w(j, i) = std::fabs(w(j, i) + (2.0 * distr(eng)) * (w_target(j,0) - w(j, i)));
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
					x_pattern(i, newDistr(eng)) = LB + (UP - LB) * distr(eng);
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
						x_pattern(i, p) = x_pattern(i, p) + (2.0*distr(eng))*(x_target(0, p) - x_pattern(i, p));
						//equation 13, page 751
					}
			}
			//return x_pattern, which is satisfy the restriction. stop 
		}
		//step 10: calculate the objective function 
		for (size_t q = 0; q < NPOP; q++)
		{
			cost(q, 0) = f(x_pattern.row(q));
		}
		//return all cost. stop 

		//step 11: update the target solution 
		//find minimum cost
		for (size_t r = 0; r < NPOP; r++)
		{
			if (cost(r,0)==cost.minCoeff()) {
				value_index1(0, 0) = cost(r, 0);
				value_index1(0, 1) = r;
				break;
			}
		}

		//find max cost
		for (size_t s = 0; s < NPOP; s++) 
		{
			if (cost(s, 0) == cost.maxCoeff()) {
				value_index2(0, 0) = cost(s, 0);
				value_index2(0, 1) = s;
				break;
			}
		}
		if (value_index1(0,0)<target )
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
		std::cout << target << std::endl;

		//return new beta. stop 
		++begin;
	}
	std::cout << target << std::endl;
	std::cout << "xtarge \t: " << x_target << std::endl;
}
double f(Eigen::Matrix<double, 1, nvars> x ){
	double PENALTY = std::pow(10, 15.0);
	int caseNum =12;
	double c[NPOP];
	double sumConstrains = 0.0;
	double temp = 0;
	double p = 1.0;
	int n = 2;

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
			temp += std::sin(x(0, i)) * std::pow( (std::sin(i+1 * std::pow(x(0, i), 2.0)/M_PI)), 2.0*10.0);
		}
		return -temp;
		break;
	case 4:
		//sum square
		temp = 0.0;
		for (size_t i = 0; i < nvars; i++)
		{
			temp = temp+(i + 1) * std::pow(x(0, i), 2.0);
		}
		return temp;
		break;
	case 5:
		//rosenbrock function
		temp = 0.0;
		for (size_t i = 1; i < nvars; i++)
		{
			temp = temp + 100.0 * (std::pow(x(0, i ) - std::pow(x(0, i-1), 2.0), 2.0)) +std::pow(1.0-x(0,i-1),2.0);
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
		return 10.0*nvars+temp;
		break;
	case 7:
		//hump function
		temp = 1.0316285 + 4 * std::pow(x(0,0),2.0) - 
			2.1 * std::pow(x(0,0),4.0) + std::pow(x(0,0),6.0) / 3 + x(0,0) * x(0,1) - 4 * std::pow(x(0,1), 2.0) + 4 * std::pow(x(0,1),4.0);
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
		return (temp + PENALTY * std::pow(std::max(sumConstrains,0.0),2.0));
		break ;
	case 10:
		//page 417. stewart calculus, question no 19 : 0.707822 0.353192
		temp = std::exp(-x(0, 0) * x(0, 1));
		sumConstrains = std::pow(x(0, 0), 2.0) + 4.0*std::pow(x(0, 1), 2.0) - 1.0;
		return (temp + PENALTY  * std::pow(std::max(sumConstrains, 0.0), 2.0));
		break;
	case 11:
		//simple square case : 1.95418 0.0885032, page 425
		temp = (1 / 2.0) * std::pow(x(0, 0) - 2.0, 2.0) + (1/2.0)*std::pow(x(0, 1) - 1 / 2.0, 2.0);
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
	default:
		break;
	}
	return 0;
}
