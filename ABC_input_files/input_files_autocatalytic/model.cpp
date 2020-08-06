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
	const double K_1 = part_params[0];
	const double K_2 = part_params[1];
	const double K_3 = part_params[2];
	const double K_4 = part_params[3];
	const double m = part_params[4];

	//Species order is: x_1 x_2 x_3 x_4 
	dydt[0] = y[0]*(-K_1*y[1]*(m - y[1]) + K_4*y[3]);
	dydt[1] = y[1]*(K_1*y[0]*(m - y[1]) - K_2*y[2]);
	dydt[2] = y[2]*(K_2*y[1] - K_3*y[3]);
	dydt[3] = y[3]*(K_3*y[2] - K_4*y[0]);

}

void Models::jac_0(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double K_1 = part_params[0];
	const double K_2 = part_params[1];
	const double K_3 = part_params[2];
	const double K_4 = part_params[3];
	const double m = part_params[4];

	J( 0 , 0 ) = -K_1*y[1]*(m - y[1]) + K_4*y[3];
	J( 0 , 1 ) = y[0]*(K_1*y[1] - K_1*(m - y[1]));
	J( 0 , 2 ) = 0;
	J( 0 , 3 ) = K_4*y[0];
	J( 1 , 0 ) = K_1*y[1]*(m - y[1]);
	J( 1 , 1 ) = -K_1*y[0]*y[1] + K_1*y[0]*(m - y[1]) - K_2*y[2];
	J( 1 , 2 ) = -K_2*y[1];
	J( 1 , 3 ) = 0;
	J( 2 , 0 ) = 0;
	J( 2 , 1 ) = K_2*y[2];
	J( 2 , 2 ) = K_2*y[1] - K_3*y[3];
	J( 2 , 3 ) = -K_3*y[2];
	J( 3 , 0 ) = -K_4*y[3];
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = K_3*y[3];
	J( 3 , 3 ) = K_3*y[2] - K_4*y[0];

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;

}
