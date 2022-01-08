#include <iostream>
#include <Eigen/dense>
#include <random> 
#include <algorithm>
#define NPOP 100
#define MAX  100
#define LB 0 
#define UP 1.0
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

	while (begin!=1)
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

		//step 12: update beta, equation 11
		beta *= 0.99;
		if (beta<0.01)
		{
			beta = 0.05;
		}
		//return new beta. stop 
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
		++begin;
	}
	std::cout << "xtarge \t: " << x_target << std::endl;
}
double f(Eigen::Matrix<double, 1, nvars> x ){
	return (x(0, 0) + x(0, 1));
}
