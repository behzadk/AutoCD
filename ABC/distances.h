// Header guard
#ifndef __DISTANCE_H_INCLUDED__
#define __DISTANCE_H_INCLUDED__
#include "particle_sim_opemp.h"


class DistanceFunctions {
	
	public:
		DistanceFunctions();
		std::vector<double> extract_species_to_fit(std::vector<state_type>&, int, int);
		long double calculate_final_gradient(std::vector<double>&);
		std::vector<double> get_signal_gradient(std::vector<double> &);
		double standard_deviation(std::vector<double>&);
		bool has_negative_species(std::vector<state_type>&);
		
		void test_fft(float, int, int, float);
		void fft_freq(double[], int, float);

		std::vector<double> arange(int, int, float);

		double get_period_frequency(std::vector<double>&, const float);
		void find_signal_peaks_and_troughs(std::vector<double>&, std::vector<int>&, std::vector<int>&);
		std::vector<double> get_amplitudes(std::vector<double>&, std::vector<int>&, std::vector<int>&);

		double get_sum_stdev(std::vector<state_type>&, int, int);
		long double get_sum_grad(std::vector<state_type>&, int);
		boost::python::list get_all_species_grads(std::vector<state_type>&, int);
		std::vector<std::vector<double>> stable_dist(std::vector<state_type>&, std::vector<int>, bool);
		std::vector<std::vector<double>> osc_dist(std::vector<state_type>&, std::vector<int>, bool, const float dt);
		std::vector<std::vector<double>> survival_dist(std::vector<state_type>&, std::vector<int>, bool);

};


#endif
