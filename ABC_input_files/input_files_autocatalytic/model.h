// Header guard
#ifndef __MODELS_H_INCLUDED__
#define __MODELS_H_INCLUDED__
#include <array>
#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include "model.h"

class Models{
	typedef boost::numeric::ublas::vector< double > ublas_vec_t;
	typedef boost::numeric::ublas::matrix< double > ublas_mat_t;
	typedef void (Models::*model_t)(const std::vector<double> &, std::vector<double> &, double, std::vector<double>&);
	typedef void (Models::*model_ublas_t)(const ublas_vec_t  &, ublas_vec_t &, double, std::vector<double>&);
	typedef void (Models::*model_jac_t)(const ublas_vec_t & x , ublas_mat_t &J , const double &, ublas_vec_t &dfdt,  std::vector<double>&);

	public:
		Models();
		std::vector < model_t > models_vec;
		std::vector < model_ublas_t > models_ublas_vec;

		std::vector < model_jac_t > models_jac_vec;
		void run_model_ublas(const ublas_vec_t  &, ublas_vec_t &, double, std::vector <double>&, int&);
		void run_jac(const ublas_vec_t &, ublas_mat_t &, const double & , ublas_vec_t &, std::vector <double> &, int &);

		void model_0(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_0(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

};

#endif