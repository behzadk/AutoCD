// Header guard
#ifndef __PARTICLE_SIM_OPEMP_H_INCLUDED__
#define __PARTICLE_SIM_OPEMP_H_INCLUDED__

// Forward declared dependencies (none)
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/python.hpp>
#include <boost/python/args.hpp>
#include <boost/thread/thread.hpp>
#include <omp.h>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include "model.h"

// Include dependencies
using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::python;


// The class
typedef std::vector<double> state_type;
typedef std::vector<double> time_type;
typedef runge_kutta_dopri5<
          state_type , double ,
          state_type , double ,
          openmp_range_algebra
        > dopri5_stepper_type;


// typedef rosenbrock4<double, double, double> stepper_type;
typedef void (Models::*model_t)(const std::vector<double> &, std::vector<double> &, double, std::vector<double>&);

typedef boost::numeric::ublas::vector< double > ublas_vec_t;
typedef boost::numeric::ublas::matrix< double > ublas_mat_t;



class Particle
{
	private:
		Models m;
		std::vector<std::vector<double>> sim_distances;

	public:
		// Constructor
		Particle(ublas_vec_t, std::vector <double>, Models, int);

		model_t particle_model;
		int particle_model_int;

		int model_ref;

		std::vector <double> part_params;
		ublas_vec_t state_init;
		vector <state_type> state_vec;
		time_type times_array;

		void set_distance_vector(std::vector<std::vector<double>>);
		std::vector<std::vector<double>> get_sim_distances() {return sim_distances;};


		void simulate_particle_rosenbrock(std::vector<double>, double, double, std::vector<int>);

		// Overloaded operators to run model step and jacobian
		void operator() ( const ublas_vec_t & , ublas_vec_t &, double); //  model
		void operator() ( const ublas_vec_t & , ublas_mat_t & ,const double &, ublas_vec_t &); // jac

		std::vector<state_type>& get_state_vec();
		boost::python::list get_state_pylist();
		boost::python::list get_end_state_jacobian();
		boost::python::list get_init_state_jacobian();
		boost::python::list get_jacobian(boost::python::list);

		boost::python::list get_final_species_values();
		double get_sum_stdev(int);
		// long double get_sum_grad();
		boost::python::list get_all_grads();


		boost::python::list py_model_func(boost::python::list);


		double get_trace();
		
		bool integration_failed = false;
		std::string integration_error;
};


// Converts a vector to boost::python::list
template <class T>
void vec_to_pylist(std::vector<T> vector) {
    typename std::vector<T>::iterator iter;
    boost::python::list list;
    for (iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(iter);
    }
}


#endif
