#include <array>
#include <iostream>
#include <vector>
#include <math.h>
#include "model.h"

Models::Models() {
	 models_ublas_vec = {&Models::model_0, &Models::model_1, &Models::model_2};
	 models_jac_vec = {&Models::jac_0, &Models::jac_1, &Models::jac_2};
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
	const double C_extra = part_params[0];
	const double C_OD = part_params[1];
	const double D = part_params[2];
	const double g_1 = part_params[3];
	const double K_A_B_1 = part_params[4];
	const double K_mu_glu = part_params[5];
	const double k_omega_B_1 = part_params[6];
	const double kA_1 = part_params[7];
	const double kB_max_1 = part_params[8];
	const double mu_max_1 = part_params[9];
	const double n_A_B_1 = part_params[10];
	const double n_omega = part_params[11];
	const double omega_max = part_params[12];
	const double S0_glu = part_params[13];

	//Species order is: N_1 S_glu B_1 A_1 
	dydt[0] = -D*y[0] + y[0]*y[1]*mu_max_1/(K_mu_glu + y[1]) - y[0]*omega_max*std::pow(y[2]*C_extra, n_omega)/(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega));
	dydt[1] = -C_OD*y[0]*y[1]*mu_max_1/(g_1*(K_mu_glu + y[1])) + D*(S0_glu - y[1]);
	dydt[2] = -y[2]*D + C_OD*y[0]*kB_max_1*std::pow(y[3]*C_extra, n_A_B_1)/(C_extra*(std::pow(K_A_B_1, n_A_B_1) + std::pow(y[3]*C_extra, n_A_B_1)));
	dydt[3] = -y[3]*D + C_OD*y[0]*kA_1/C_extra;

}
void Models::jac_0(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double C_extra = part_params[0];
	const double C_OD = part_params[1];
	const double D = part_params[2];
	const double g_1 = part_params[3];
	const double K_A_B_1 = part_params[4];
	const double K_mu_glu = part_params[5];
	const double k_omega_B_1 = part_params[6];
	const double kA_1 = part_params[7];
	const double kB_max_1 = part_params[8];
	const double mu_max_1 = part_params[9];
	const double n_A_B_1 = part_params[10];
	const double n_omega = part_params[11];
	const double omega_max = part_params[12];
	const double S0_glu = part_params[13];

	J( 0 , 0 ) = -D + y[1]*mu_max_1/(K_mu_glu + y[1]) - omega_max*std::pow(y[2]*C_extra, n_omega)/(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega));
	J( 0 , 1 ) = -y[0]*y[1]*mu_max_1/std::pow(K_mu_glu + y[1], 2) + y[0]*mu_max_1/(K_mu_glu + y[1]);
	J( 0 , 2 ) = y[0]*n_omega*omega_max*std::pow(y[2]*C_extra, 2*n_omega)/(y[2]*std::pow(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega), 2)) - y[0]*n_omega*omega_max*std::pow(y[2]*C_extra, n_omega)/(y[2]*(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega)));
	J( 0 , 3 ) = 0;
	J( 1 , 0 ) = -C_OD*y[1]*mu_max_1/(g_1*(K_mu_glu + y[1]));
	J( 1 , 1 ) = C_OD*y[0]*y[1]*mu_max_1/(g_1*std::pow(K_mu_glu + y[1], 2)) - C_OD*y[0]*mu_max_1/(g_1*(K_mu_glu + y[1])) - D;
	J( 1 , 2 ) = 0;
	J( 1 , 3 ) = 0;
	J( 2 , 0 ) = C_OD*kB_max_1*std::pow(y[3]*C_extra, n_A_B_1)/(C_extra*(std::pow(K_A_B_1, n_A_B_1) + std::pow(y[3]*C_extra, n_A_B_1)));
	J( 2 , 1 ) = 0;
	J( 2 , 2 ) = -D;
	J( 2 , 3 ) = -C_OD*y[0]*kB_max_1*n_A_B_1*std::pow(y[3]*C_extra, 2*n_A_B_1)/(y[3]*C_extra*std::pow(std::pow(K_A_B_1, n_A_B_1) + std::pow(y[3]*C_extra, n_A_B_1), 2)) + C_OD*y[0]*kB_max_1*n_A_B_1*std::pow(y[3]*C_extra, n_A_B_1)/(y[3]*C_extra*(std::pow(K_A_B_1, n_A_B_1) + std::pow(y[3]*C_extra, n_A_B_1)));
	J( 3 , 0 ) = C_OD*kA_1/C_extra;
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;

}
void Models::model_1(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double C_extra = part_params[0];
	const double C_OD = part_params[1];
	const double D = part_params[2];
	const double g_1 = part_params[3];
	const double K_A_B_1 = part_params[4];
	const double K_mu_glu = part_params[5];
	const double k_omega_B_1 = part_params[6];
	const double kA_1 = part_params[7];
	const double kB_max_1 = part_params[8];
	const double mu_max_1 = part_params[9];
	const double n_A_B_1 = part_params[10];
	const double n_omega = part_params[11];
	const double omega_max = part_params[12];
	const double S0_glu = part_params[13];

	//Species order is: N_1 S_glu B_1 A_1 
	dydt[0] = -D*y[0] + y[0]*y[1]*mu_max_1/(K_mu_glu + y[1]) - y[0]*omega_max*std::pow(y[2]*C_extra, n_omega)/(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega));
	dydt[1] = -C_OD*y[0]*y[1]*mu_max_1/(g_1*(K_mu_glu + y[1])) + D*(S0_glu - y[1]);
	dydt[2] = -y[2]*D + C_OD*std::pow(K_A_B_1, n_A_B_1)*y[0]*kB_max_1/(C_extra*(std::pow(K_A_B_1, n_A_B_1) + std::pow(y[3]*C_extra, n_A_B_1)));
	dydt[3] = -y[3]*D + C_OD*y[0]*kA_1/C_extra;

}
void Models::jac_1(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double C_extra = part_params[0];
	const double C_OD = part_params[1];
	const double D = part_params[2];
	const double g_1 = part_params[3];
	const double K_A_B_1 = part_params[4];
	const double K_mu_glu = part_params[5];
	const double k_omega_B_1 = part_params[6];
	const double kA_1 = part_params[7];
	const double kB_max_1 = part_params[8];
	const double mu_max_1 = part_params[9];
	const double n_A_B_1 = part_params[10];
	const double n_omega = part_params[11];
	const double omega_max = part_params[12];
	const double S0_glu = part_params[13];

	J( 0 , 0 ) = -D + y[1]*mu_max_1/(K_mu_glu + y[1]) - omega_max*std::pow(y[2]*C_extra, n_omega)/(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega));
	J( 0 , 1 ) = -y[0]*y[1]*mu_max_1/std::pow(K_mu_glu + y[1], 2) + y[0]*mu_max_1/(K_mu_glu + y[1]);
	J( 0 , 2 ) = y[0]*n_omega*omega_max*std::pow(y[2]*C_extra, 2*n_omega)/(y[2]*std::pow(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega), 2)) - y[0]*n_omega*omega_max*std::pow(y[2]*C_extra, n_omega)/(y[2]*(std::pow(k_omega_B_1, n_omega) + std::pow(y[2]*C_extra, n_omega)));
	J( 0 , 3 ) = 0;
	J( 1 , 0 ) = -C_OD*y[1]*mu_max_1/(g_1*(K_mu_glu + y[1]));
	J( 1 , 1 ) = C_OD*y[0]*y[1]*mu_max_1/(g_1*std::pow(K_mu_glu + y[1], 2)) - C_OD*y[0]*mu_max_1/(g_1*(K_mu_glu + y[1])) - D;
	J( 1 , 2 ) = 0;
	J( 1 , 3 ) = 0;
	J( 2 , 0 ) = C_OD*std::pow(K_A_B_1, n_A_B_1)*kB_max_1/(C_extra*(std::pow(K_A_B_1, n_A_B_1) + std::pow(y[3]*C_extra, n_A_B_1)));
	J( 2 , 1 ) = 0;
	J( 2 , 2 ) = -D;
	J( 2 , 3 ) = -C_OD*std::pow(K_A_B_1, n_A_B_1)*y[0]*kB_max_1*n_A_B_1*std::pow(y[3]*C_extra, n_A_B_1)/(y[3]*C_extra*std::pow(std::pow(K_A_B_1, n_A_B_1) + std::pow(y[3]*C_extra, n_A_B_1), 2));
	J( 3 , 0 ) = C_OD*kA_1/C_extra;
	J( 3 , 1 ) = 0;
	J( 3 , 2 ) = 0;
	J( 3 , 3 ) = -D;

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	dfdt[2] = 0.0;
	dfdt[3] = 0.0;

}
void Models::model_2(const ublas_vec_t  &y , ublas_vec_t &dydt , double t, std::vector <double> &part_params)
{
	//Unpack parameters
	const double C_OD = part_params[0];
	const double D = part_params[1];
	const double g_1 = part_params[2];
	const double K_mu_glu = part_params[3];
	const double mu_max_1 = part_params[4];
	const double S0_glu = part_params[5];

	//Species order is: N_1 S_glu 
	dydt[0] = -D*y[0] + y[0]*y[1]*mu_max_1/(K_mu_glu + y[1]);
	dydt[1] = -C_OD*y[0]*y[1]*mu_max_1/(g_1*(K_mu_glu + y[1])) + D*(S0_glu - y[1]);

}
void Models::jac_2(const ublas_vec_t & y , ublas_mat_t &J , const double &/* t*/ , ublas_vec_t &dfdt, std::vector <double> &part_params)
{
	//Unpack parameters
	const double C_OD = part_params[0];
	const double D = part_params[1];
	const double g_1 = part_params[2];
	const double K_mu_glu = part_params[3];
	const double mu_max_1 = part_params[4];
	const double S0_glu = part_params[5];

	J( 0 , 0 ) = -D + y[1]*mu_max_1/(K_mu_glu + y[1]);
	J( 0 , 1 ) = -y[0]*y[1]*mu_max_1/std::pow(K_mu_glu + y[1], 2) + y[0]*mu_max_1/(K_mu_glu + y[1]);
	J( 1 , 0 ) = -C_OD*y[1]*mu_max_1/(g_1*(K_mu_glu + y[1]));
	J( 1 , 1 ) = C_OD*y[0]*y[1]*mu_max_1/(g_1*std::pow(K_mu_glu + y[1], 2)) - C_OD*y[0]*mu_max_1/(g_1*(K_mu_glu + y[1])) - D;

	dfdt[0] = 0.0;
	dfdt[1] = 0.0;

}
