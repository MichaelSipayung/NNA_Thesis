#include <iostream>
#include <Eigen/dense>
#include <random> 
#define NPOP 100
#define MAX  100
#define LB 0 
#define UP 1.0
#define nvars 2
double f(Eigen::Matrix<double, 1, nvars>);

int main()
{
	Eigen::Matrix<double, NPOP, 1> x_lb;
	Eigen::Matrix<double, NPOP, 1> x_ub;
	Eigen::Matrix<double, NPOP, 1> cost;
	Eigen::Matrix<double, NPOP, nvars> x_pattern;
	Eigen::Matrix<double, 1, NPOP> ww; 
	Eigen::Matrix<double, NPOP, NPOP> w;
	Eigen::Matrix<double, 1, NPOP-1> t;
	Eigen::Matrix <double, 1, 2> value_index; 
	double beta = 1.0; 
	constexpr int FLOAT_MIN = 0;
	constexpr int FLOAT_MAX = 1;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);
	
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
	// 
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
			if (w(var, j) == 0.0) {
				w(var, j) = t(ink);
				++ink;
			}
		}
	}
	//step 4: return weight matrix, which  is satisfy the constrains. stop 
}

double f(Eigen::Matrix<double, 1, nvars> x ){
	return (x(0, 0) + x(0, 1));
}
