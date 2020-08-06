// Header guard
#ifndef __POPULATION_H_INCLUDED__
#define __POPULATION_H_INCLUDED__
#include "particle_sim_opemp.h"


typedef std::vector<std::vector<std::vector<double>>> PopDistances;
class Population {
	public:
		Population(const int, const float, const float, 
			const float, boost::python::list, boost::python::list, boost::python::list, boost::python::list, double, double);

		std::vector< std::vector<double> > unpack_parameters(boost::python::list);
		std::vector<int> unpack_model_references(boost::python::list);
		std::vector<int> fit_species;

		void generate_particles();
		void simulate_particles();
		void calculate_particle_distances(int);
		void accumulate_distances();
		
		std::vector< ublas_vec_t > unpack_parameters_to_ublas(boost::python::list);

		PopDistances get_population_distances() {return _all_distances;};
		
		boost::python::list get_flattened_distances_list();
		boost::python::list get_timepoints_list();
		boost::python::list get_particle_state_list(int);
		// boost::python::list get_particle_eigenvalues(int);
		boost::python::list get_particle_init_state_jacobian(int);
		boost::python::list get_particle_end_state_jacobian(int);
		boost::python::list get_particle_final_species_values(int);
		boost::python::list get_particle_jacobian(boost::python::list, int);

		boost::python::list py_model_func(boost::python::list, int);
		double get_particle_sum_stdev(int, int);
		// long double get_particle_sum_grad(int);
		boost::python::list get_particle_grads(int);

		// double get_particle_det(int);
		// void get_particle_laplace_expansion(int);

		bool check_integration_failure(int);
		std::string get_particle_integration_error(int);
		boost::python::list get_all_particle_integration_errors();
		double get_particle_trace(int);

	private:
		int _n_sims;
		float  _t_0;
		float _t_end;
		float _dt;
		int _distance_function_mode;
		double _abs_tol;
		double _rel_tol;

		std::vector<double> _time_array;
		std::vector<int> _model_refs;
		PopDistances _all_distances;
		std::vector<std::vector<double>> _all_params;
		std::vector<Particle> _particle_vector;
		std::vector<ublas_vec_t>_all_state_init;

};




#endif
