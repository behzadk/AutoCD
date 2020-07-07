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

		void model_1(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_1(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_2(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_2(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_3(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_3(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_4(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_4(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_5(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_5(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_6(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_6(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_7(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_7(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_8(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_8(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_9(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_9(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_10(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_10(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_11(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_11(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_12(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_12(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_13(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_13(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_14(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_14(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_15(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_15(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_16(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_16(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_17(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_17(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_18(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_18(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_19(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_19(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_20(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_20(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_21(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_21(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_22(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_22(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_23(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_23(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_24(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_24(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_25(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_25(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_26(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_26(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_27(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_27(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_28(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_28(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_29(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_29(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_30(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_30(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_31(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_31(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_32(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_32(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_33(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_33(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_34(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_34(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_35(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_35(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_36(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_36(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_37(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_37(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_38(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_38(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_39(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_39(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_40(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_40(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_41(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_41(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_42(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_42(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_43(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_43(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_44(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_44(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_45(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_45(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_46(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_46(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_47(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_47(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_48(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_48(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_49(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_49(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_50(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_50(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_51(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_51(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_52(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_52(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_53(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_53(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_54(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_54(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_55(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_55(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_56(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_56(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_57(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_57(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_58(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_58(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_59(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_59(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_60(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_60(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_61(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_61(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_62(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_62(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_63(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_63(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_64(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_64(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_65(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_65(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_66(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_66(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_67(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_67(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

		void model_68(const ublas_vec_t  & , ublas_vec_t & , double , std::vector <double> &);
		void jac_68(const ublas_vec_t &, ublas_mat_t & , const double &, ublas_vec_t &, std::vector <double> &);

};

#endif