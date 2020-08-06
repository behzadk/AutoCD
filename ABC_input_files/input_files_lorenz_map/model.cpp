#include <array>
#include <iostream>
#include <vector>
#include <math.h>
#include "model.h"

Models::Models() {
	 models_ublas_vec = {&Models::model_0};
	 models_jac_vec = {&Models::jac_0};
};

void Models::run_model_ublas(const ublas_vec_t &y , ublas_vec_t &dxdt , double t, std::vector <double> &part_params, int &model_ref)
{
	(this->*models_ublas_vec[model_ref])(y, dxdt, t, part_params);
}

void Models::run_jac(const ublas_vec_t & x , ublas_mat_t &J , const double & t , ublas_vec_t &dfdt, std::vector <double> &part_params, int &model_ref)
{
	(this->*models_jac_vec[model_ref])(x, J, t, dfdt, part_params);
}

void Models::model_0(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double beta_1 = part_params[0];
	const double p = part_params[1];
	const double sigma = part_params[2];

	//Species order is: X_1 Y_1 Z_1 
	dydt[0] = sigma*(-y[0] + y[1]);
	dydt[1] = y[0]*(-y[2] + p) - y[1];
	dydt[2] = y[0]*y[1] - y[2]*beta_1;

}
void Models::jac_0(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double beta_1 = part_params[0];
	const double p = part_params[1];
	const double sigma = part_params[2];

	J( 0 , 0 ) = -sigma;
	J( 0 , 1 ) = sigma;
	J( 0 , 2 ) = 0;
	J( 1 , 0 ) = -y[2] + p;
	J( 1 , 1 ) = -1;
	J( 1 , 2 ) = -y[0];
	J( 2 , 0 ) = y[1];
	J( 2 , 1 ) = y[0];
	J( 2 , 2 ) = -beta_1;

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;

}
