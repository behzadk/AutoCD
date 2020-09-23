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
	const double alpha_1 = part_params[0];
	const double alpha_2 = part_params[1];
	const double beta_1 = part_params[2];
	const double beta_2 = part_params[3];
	const double d = part_params[4];
	const double K_1 = part_params[5];
	const double K_2 = part_params[6];
	const double rho_1 = part_params[7];
	const double rho_2 = part_params[8];

	//Species order is: N_1 N_2 P_1 
	dydt[0] = -y[0]*y[2]*alpha_1 + y[0]*rho_1*(1 - y[0]/K_1);
	dydt[1] = -y[1]*y[2]*alpha_2 + y[1]*rho_2*(1 - y[1]/K_2);
	dydt[2] = y[2]*(y[0]*beta_1 + y[1]*beta_2 - d);

}
void Models::jac_0(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double alpha_1 = part_params[0];
	const double alpha_2 = part_params[1];
	const double beta_1 = part_params[2];
	const double beta_2 = part_params[3];
	const double d = part_params[4];
	const double K_1 = part_params[5];
	const double K_2 = part_params[6];
	const double rho_1 = part_params[7];
	const double rho_2 = part_params[8];

	J( 0 , 0 ) = -y[2]*alpha_1 + rho_1*(1 - y[0]/K_1) - y[0]*rho_1/K_1;
	J( 0 , 1 ) = 0;
	J( 0 , 2 ) = -y[0]*alpha_1;
	J( 1 , 0 ) = 0;
	J( 1 , 1 ) = -y[2]*alpha_2 + rho_2*(1 - y[1]/K_2) - y[1]*rho_2/K_2;
	J( 1 , 2 ) = -y[1]*alpha_2;
	J( 2 , 0 ) = y[2]*beta_1;
	J( 2 , 1 ) = y[2]*beta_2;
	J( 2 , 2 ) = y[0]*beta_1 + y[1]*beta_2 - d;

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;

}
