#include <array>
#include <iostream>
#include <vector>
#include <math.h>
#include "distances.h"
#include <algorithm>
#include <numeric>

extern "C" {
	#include "kissfft/kiss_fft.h"
	#include "kissfft/kiss_fftr.h"
}

DistanceFunctions::DistanceFunctions() {
}

std::vector<double> DistanceFunctions::arange(int start, int stop, float step) {
    std::vector<double> values;

    for (double value = start; value < stop; value += step) {
        values.push_back(value);
    }

    return values;
}


/*! \brief Return the Discrete Fourier Transform sample frequencies.
 *         
 *
 *  
 *	See np.fft.fftfreq
 *	
 */
void DistanceFunctions::fft_freq(double * results, int n, float step_size) {
    float val = 1.0 / (n * step_size);
    // double results[n];

    int N;
    if (n % 2 == 0) {
	    N = (n) / 2 - 1 + 1;
    }
    else {
	    N = (n-1)/2 + 1;
    }

    std::vector<double> p1 = arange(0, N, 1);
    
    for (int i=0; i <= N; i++) {
	    results[i] = p1[i] * val;
    }

    std::vector<double> p2 = arange(-N, 0, 1);
	for (int i=0; i < p2.size(); i++) {
	    results[N+i] = p2[i] * val;
    }
}

/*! \brief Returns the period frequency of a signal
 *  \param signal vector of the signal to be anaylsed
 *	\param dt float of the time step used in the simulation
 */
double DistanceFunctions::get_period_frequency(std::vector<double>& signal, const float dt)
{
	int n_samples = signal.size();

  	kiss_fft_scalar signal_scalar[n_samples];
	for (int i = 0; i < n_samples; i++) {
		signal_scalar[i] = signal[i];
	}

  	kiss_fft_cpx out[(n_samples /2) + 1];
    double period;

  	kiss_fftr_cfg cfg;

  	if ((cfg = kiss_fftr_alloc(n_samples, 0/*is_inverse_fft*/, NULL, NULL)) != NULL) 
  	{
	    size_t x;

	    kiss_fftr(cfg, signal_scalar, out);
	    free(cfg);

	    double max_real_part = pow((out[1].r + out[1].i), 2);
	   	int max_arg = 1;

	    for (x = 2; x < (n_samples/2 ) + 1; x++)
	    {
	    	double mag  = pow((out[x].r + out[x].i), 2);
	    	if (mag > max_real_part) {
	    		if ( std::isinf(mag)) {
	    			continue;
	    		}
	    		else {
		    		max_real_part = mag;
				   	max_arg = x;
		    	}
	    	}
	    }

	    double freq_bins[n_samples];
	    fft_freq(freq_bins, n_samples, dt);

	   	double max_freq = 2 * freq_bins[max_arg];
	   	period = 1/max_freq;

   	}

    return period;
}

/*! Gets the signal of a single species from a simulation.
 *  \param state_vec vector containing simulation states.
 *	\param species_idx index of species to extract.
 *  \return species_val_vec vector of simulation signal for a single species.
 */
std::vector<double> DistanceFunctions::extract_species_to_fit(std::vector<state_type>& state_vec, int species_idx, int from_time_index=0)
{
	std::vector<double> species_val_vec;

	for (auto tp_iter = state_vec.begin() + from_time_index; tp_iter != state_vec.end(); tp_iter++) {
        state_type sim_vec = *tp_iter;
        species_val_vec.push_back(sim_vec[species_idx]);
    }

    return species_val_vec;
}

/*! Generates a derivative of a signal vector.
 *  \param signal vector states for a single species
 *  \return signal_gradient vector containing gradients at each time point
 */
std::vector<double> DistanceFunctions::get_signal_gradient(std::vector<double>& signal)
{
	std::vector<double> signal_gradient;

	for (int i = 1; i < signal.size(); i++) {
		double grad = signal[i] - signal[i-1];
		signal_gradient.push_back(grad);
	}

	return signal_gradient;

}
/*! Gets the gradient at the final timepoint of a species simulation.
 *  \param species_vals vector containing the states of a single species
 *  \return grad the gradient at the final time point
 */
long double DistanceFunctions::calculate_final_gradient(std::vector<double>& species_vals)
{
	double i = species_vals.end()[-2];
	double j = species_vals.end()[-1];

	long double grad = j - i;
	return grad;
}


/*!
 *	\brief Appends the indexes of each peak and trough to the supplied vectors
 *  Iterates through the signal gradient identifying changes between positive and negative
 *	gradients. Differences in sign of the two current and previous gradient indicate a
 *	peak or a trough.
 *
 *  \param signal_gradient vector containing the gradients of a simulation
 *  \param peak_idx vector in which the identified peak indexes are stored
 *  \param trough_idx vector in which the identified trough indexes are stored
 *
 */
void DistanceFunctions::find_signal_peaks_and_troughs(std::vector<double>& signal_gradient, std::vector<int>& peak_idx, std::vector<int>& trough_idx)
{
	double current_gradient;
	double previous_gradient;
	
	for (int i = 1; i < signal_gradient.size(); i++) {
		// Set current and previous gradients
		current_gradient = signal_gradient[i];
		previous_gradient = signal_gradient[i-1];

		/* 
		A peak is when positive state precedes negative state.
		A trough is when negative state precedes positive state.
		*/
		if (current_gradient < 0 && previous_gradient > 0) {
			peak_idx.push_back(i);
		} else if (current_gradient > 0 && previous_gradient < 0) {
			trough_idx.push_back(i);
		} else {
			continue;
		}
	}
}

/*! \brief Calculates amplitudes from signal and precalculated peaks and trough indexes.
 *
 *  Iterates through the precomputed peak and trough indexes, calculating the amplitudes between them
 *  \param signal vector containing the states of a simulated species
 *  \param peak_idx vector containing indexes of peaks
 *  \param trough_idx vector containing indexes of troughs
 *
 *  \return amplitudes vector of amplitudes in a signal.
 */
std::vector<double> DistanceFunctions::get_amplitudes(std::vector<double>& signal, std::vector<int>& peak_idx, std::vector<int>& trough_idx) 
{
	int num_peaks = peak_idx.size();
	int num_troughs = trough_idx.size();

	int max_iterations = min(num_peaks, num_troughs);

	std::vector<double> amplitudes;

	for (int i = 0; i < max_iterations; i++)
	{	
		int peak_i = peak_idx[i];
		int trough_i = trough_idx[i];

		double amp = abs(signal[peak_i] - signal[trough_i]);
		amplitudes.push_back(amp);
	}

	return amplitudes;
}

/*! \brief Calculates standard deviation of a signal.
 *        
 *	Calculates standard deviation of a signal.
 *  \param signal vector containing signal for which std is calculated
 *  \return stdev the standard deviation
 */	
double DistanceFunctions::standard_deviation(std::vector<double>& signal) {
	double mean;
	double sum_signal = 0.0;

	for (int i = 0; i < signal.size(); i++) {
		double val = signal[i];
		sum_signal = sum_signal + signal[i];
	}

	mean = sum_signal / signal.size();

	double sq_sum = 0.0;

	//sum square
	for (int i = 0; i < signal.size(); i++) {
		double val = std::pow( (signal[i] - mean) , 2);
		sq_sum += val;
	}

	if (sq_sum == 0){
		return 0;
	}
	double stdev = sqrt( sq_sum / (signal.size()) );

	return stdev;
}

/*! \brief Returns bool indicating of the simulation contains negative species.
 *  Iterates through all indexes of a vector, with each index containing n-datapoints (one for each species).
 *  Gets the minimum value at each timepoint, returning True if the minimum value is less than 0.
 *
 *	\param state_vec vector containing the states of all species for each timepoint.
 *  \return bool True if the simulation contains negative species.
 */
bool DistanceFunctions::has_negative_species(std::vector<state_type>& state_vec) {
	for (auto tp_iter = state_vec.begin(); tp_iter != state_vec.end(); tp_iter++) {
		std::vector<double> sim_vec = *tp_iter;
		double min_value = *std::min_element(sim_vec.begin(),sim_vec.end());

		if (min_value < 0) {
			return true;
		}
	}
	return false;
}

/*! \brief Calculates the sum of standard deviations for all species.
 *  Iterates through the species to fit and caluclates the standard deviation for each. The sum of the
 *  standard deviations are returned.
 *	\param state_vec vector containing states of all species for each timepoint in a simulation
 *  \param n_species integer indicating the number of species in the simulation.
 *  \return sum_stdev The sum of standard deviations for species in the simulation.
 */
double DistanceFunctions::get_sum_stdev(std::vector<state_type>& state_vec, int n_species, int from_time_index) {
	double sum_stdev = 0;

	for (int i = 0; i < n_species; i++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, i, from_time_index);
		double stdev = standard_deviation(signal);
		sum_stdev = sum_stdev + stdev;
	}

	return sum_stdev;
}

long double DistanceFunctions::get_sum_grad(std::vector<state_type>& state_vec, int n_species) {

	long double sum_grad = 0;
	for (int i = 0; i < n_species; i++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, n_species, 0);

		long double final_gradient = calculate_final_gradient(signal);
		sum_grad = sum_grad + final_gradient;
	}

	return sum_grad;

}

boost::python::list DistanceFunctions::get_all_species_grads(std::vector<state_type>& state_vec, int n_species) {
	boost::python::list all_grads;

	for (int i = 0; i < n_species; i++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, n_species, 0);

		std::vector<double> signal_gradient = get_signal_gradient(signal);

		long double final_gradient = calculate_final_gradient(signal);
		all_grads.append(final_gradient);
	}
	return all_grads;
}

/*! \brief Calculates distances for oscillatory objective.
 *  WARNING: CONTAINS HARDCODED PARAMETERS
 *
 */	
std::vector<std::vector<double>> DistanceFunctions::osc_dist(std::vector<state_type>& state_vec, std::vector<int> species_to_fit, bool integration_failed, const float dt) {
	std::vector<std::vector<double>> sim_distances;

	double max_dist = std::numeric_limits<double>::max();

	std::vector<double> max_distances = {max_dist, max_dist, max_dist};

	// If integration has failed, return all distances as maximum
	if (integration_failed) {
		for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
			sim_distances.push_back(max_distances);
		}
		return sim_distances;
	}

	int from_time_index = 50;

	int amplitude_threshold = 0;

	for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, *it, from_time_index);
		std::vector<double> signal_gradient = get_signal_gradient(signal);
		
		std::vector<int> peak_indexes;
		std::vector<int> trough_indexes;

		find_signal_peaks_and_troughs(signal_gradient, peak_indexes, trough_indexes);
		
		double signal_period_freq = get_period_frequency(signal, dt);

		std::vector<double> singal_amplitudes = get_amplitudes(signal, peak_indexes, trough_indexes);

		double threshold_amplitudes_count;
		double final_amplitude;

		// If no amplitudes, set threshold and final amplitude to zero
		if (singal_amplitudes.size() == 0) {
			threshold_amplitudes_count = 0;
			final_amplitude = 0;
		}
		else {
			threshold_amplitudes_count = 0;

			// Count threshold amplitudes
			for (auto amp_it = singal_amplitudes.begin(); amp_it != singal_amplitudes.end(); amp_it++) {
				if (*amp_it > amplitude_threshold) {
					threshold_amplitudes_count = threshold_amplitudes_count + 1;
				}
			}

			// Set final amplitude
			final_amplitude = singal_amplitudes[singal_amplitudes.size() - 1];
		}

		std::vector<double> signal_distances = {threshold_amplitudes_count, final_amplitude, signal_period_freq};
		sim_distances.push_back(signal_distances);

	}

	return sim_distances;
}


/*! \brief Calculates distances for stable objective.
 *  Iterates the species that should be fit, returning a vector with the final gradient, standard deviation and
 *  final species value. Standard deviation is calculated from the final 10 percent of the simulation, this is currently
 *  hard coded.
 *
 *  \param state_vec vector containing the states of a simulation for all speices
 *  \param species_to_fit vector containing the indexes of species that should be fit.
 *  \param integration_failed bool indicating if the integration failed
 *  \param sim_distances vector in which the calculated distances should be stored

 *  \return sim_distances vector containing the calculated distances.
 */
std::vector<std::vector<double>> DistanceFunctions::stable_dist(std::vector<state_type>& state_vec, std::vector<int> species_to_fit, bool integration_failed) {
	std::vector<std::vector<double>> sim_distances;

	double max_dist = std::numeric_limits<double>::max();

	std::vector<double> max_distances = {max_dist, max_dist, max_dist};

	if (integration_failed) {
		for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
			sim_distances.push_back(max_distances);
		}

		return sim_distances;
	}

	auto sim_size = state_vec.size();
	int from_time_index = floor(sim_size - (sim_size * 0.1));

	for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, *it, from_time_index);
	}

	// Iterate through all species to fit. Extract data.
	for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, *it, from_time_index);

		std::vector<double> signal_gradient = get_signal_gradient(signal);

		double stdev = standard_deviation(signal);
		double final_gradient = fabs(signal_gradient.end()[-1]);
		double final_value = 1/signal.end()[-1];

		std::vector<double> signal_distances = {final_gradient, stdev, final_value};

		sim_distances.push_back(signal_distances);
	}

	return sim_distances;

}

/*! \brief Calculates distances for survival objective. Returns vector of distances for each species.
 *  WARNING: CONTAINS HARDCODED PARAMETERS
 */	
std::vector<std::vector<double>> DistanceFunctions::survival_dist(std::vector<state_type>& state_vec, std::vector<int> species_to_fit, bool integration_failed) {
	std::vector<std::vector<double>> sim_distances;

	double max_dist = std::numeric_limits<double>::max();

	std::vector<double> max_distances = {max_dist};

	if (integration_failed) {
		for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
			sim_distances.push_back(max_distances);
		}

		return sim_distances;
	}

	auto sim_size = state_vec.size();
	int from_time_index = floor(sim_size - (sim_size * 0.1));

	for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, *it, from_time_index);
	}

	// Iterate through all species to fit. Extract data.
	for (auto it = species_to_fit.begin(); it != species_to_fit.end(); it++) {
		std::vector<double> signal = extract_species_to_fit(state_vec, *it, from_time_index);
		double final_value = 1/signal.end()[-1];
		std::vector<double> signal_distances = {final_value};
		sim_distances.push_back(signal_distances);
	}

	return sim_distances;

}
