#include <iostream>
#define BOOST_UBLAS_TYPE_CHECK 0

#include "model.h"
#include <boost/numeric/odeint.hpp>
#include <omp.h>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include "particle_sim_opemp.h"
using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::python;
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/numeric/odeint/integrate/max_step_checker.hpp>
#include "distances.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

using namespace std;
using namespace boost::numeric::odeint;

class ScopedGILRelease
{
public:
    inline ScopedGILRelease(){
        m_thread_state = PyEval_SaveThread();
    }
    inline ~ScopedGILRelease() {
        PyEval_RestoreThread(m_thread_state);
        m_thread_state = NULL;
    }
private:
    PyThreadState * m_thread_state;
};

/*! \brief Records states when called and stops simulations early if necessary.
*   
*   The observer is called by boost at the end of each step. We have hard coded minimum species
*   values that will stop the simulation early. These are for the fitted species only. Errors are also thrown
*   if any species (not just fitted species) becomes negative.
*
*   Minimum species values should be made adjustable and provided as a parameter in future.
*/
struct simulation_observer
{
    std::vector< state_type >& m_states;
    size_t m_steps;
    rosenbrock4_controller<rosenbrock4<double> > m_stepper_controller;
    std::vector<int> m_fit_species;

    simulation_observer( std::vector< state_type > &states, rosenbrock4_controller<rosenbrock4<double> > &stepper_controller, std::vector<int> &fit_species )
    : m_states( states ), m_stepper_controller( stepper_controller ), m_fit_species(fit_species)  { }

    void operator()( ublas_vec_t x , double t )
    {
        m_steps +=1;
        std::vector<double> new_state;
        bool species_decayed = false;
        bool negative_species = false;

        for (auto it = m_fit_species.begin(); it != m_fit_species.end(); it++) {
            if(x(*it) < 1e-10) {
                species_decayed = true;
            }
        }

        // //Test
        // if(x(0) < 1e-100) {
        //     species_decayed = true;
        // }

        // //Test
        // if(x(1) < 1e-100) {
        //     species_decayed = true;
        // }


        for(int i = 0; i < x.size(); i++){
            double val = x(i);

            // if (val < 1e-200) {
            //     species_decayed = true;
            // }

            if(val < 0) {
                negative_species = true;
            };

            new_state.push_back(val);

        }
        m_states.push_back( new_state );

        if (species_decayed) throw runtime_error("species_decayed");
        if (negative_species) throw runtime_error("negative_species");

    }
};

/*! \brief Particle to be simulated.
*   
*   Particles contain all the necessary information to simulate a model  
* 
*   \param init_state   initial state of the simulation
*   \param params       vector containing parameters for the particular model
*   \param model_obj    The model to be simulated
*   \param model_idx    The model index
*/
Particle::Particle(ublas_vec_t init_state, std::vector<double> params, Models model_obj, int model_idx)
{
	Models m = model_obj;
	part_params = params;
    model_ref = model_idx;
    state_init = init_state;
}


/*!
*   \brief Run the particle's model 
*   Overloaded function run particle model.
*   
*   \param y initial species state
*   \param dxdt time step size
*   \param t initial time
*
*/
void Particle::operator() (const ublas_vec_t & y , ublas_vec_t &dxdt , double t ) // run_model_func
{
    std::vector<model_t> m_vec = m.models_vec;
    m.run_model_ublas(y, dxdt, t, part_params, model_ref);
}

/*! \brief Run the particle's Jacobian.
*   Overloaded function run particle model.
*   
*   \param x species state
*   \param J Matrix for jacobian storage vector
*   \param t time (not used)
*   \param dfdt Explicit time derivative vector
*
*/
void Particle::operator() (const ublas_vec_t & x , ublas_mat_t &J , const double & t , ublas_vec_t &dfdt ) // run_jac
{
    m.run_jac(x, J, t, dfdt, part_params, model_ref);
}

/*! \brief Simulates particle for a given vector of time points.
*
*   Simulates particle for a given vector of time points using supplied absolute and relative tolerances (abs_tol and rel_tol).
*   the fit_species vector is used for early stopping when one or more species decays. Max step check is currently hard coded.
*   
*   \param time_points  Vector of time points at which the state will be recorded
*   \param abs_tol      Absolute tolerance to use in the error stepper
*   \param rel_tol      Relative tolerance to use in the error stepper
*   \param fit_species  Species being fitted, importantly observed for early stopping criteria
*/
void Particle::simulate_particle_rosenbrock(std::vector<double> time_points, double abs_tol, double rel_tol, std::vector<int> fit_species)
{
    auto rosen_stepper = rosenbrock4_controller< rosenbrock4< double > >( abs_tol , rel_tol);

    max_step_checker mx_step = max_step_checker(1e5);

    double dt = time_points[1] - time_points[0];

    // Check if step size is negative
    if ( ( time_points[1] - time_points[0] ) < 0 ) {
        dt = dt * -1;
    }

    integrate_times(  rosen_stepper, 
        make_pair(boost::ref( *this ), boost::ref( *this )) , 
        state_init , time_points.begin(), time_points.end() , dt, simulation_observer(state_vec, rosen_stepper, fit_species), mx_step);
}

/*! \brief Python compatible model function
*
*   Exposes the model function to python. Takes an input of species values as a list and returns the dydt values as
*   a python list.
*   
*   \param input_y python list containing initial species values
*   \return output python list containing all state values
*/
boost::python::list Particle::py_model_func(boost::python::list input_y)
{
    
    int n_species = state_init.size();
    ublas_vec_t y(n_species);

    // Load python data into ublas vector
    for(int i = 0; i < n_species; i++)
    {
        double val = boost::python::extract<double>(input_y[i]);
        y(i) = val;
    }

    ublas_vec_t dydt(n_species);
    double t;

    // Simulate model
    std::vector<model_t> m_vec = m.models_vec;
    m.run_model_ublas(y, dydt, t, part_params, model_ref);

    // Convert dydt to boost python list
    boost::python::list output;
    for(int i = 0; i < n_species; i++)
    {
        output.append(dydt[i]);
    }

    return output;
}


/*! \brief Gets Jacobian trace using the initial species values.
*   \return trace value being the sum of the diagonal of the Jacobian
*/
double Particle::get_trace()
{
    int n_species = state_init.size();

    ublas_vec_t y(n_species);
    for (int i=0; i < n_species; i++) {
        y(i) = state_vec.back()[i];
    }
    
    // Init matrix n_species x n_species
    ublas_mat_t J (n_species, n_species);

    // Dummy parameters
    const double t = 0;

    ublas_vec_t dfdt(n_species);

    // Fill jacobian matrix
    m.run_jac(y, J, t, dfdt, part_params, model_ref);

    double trace = 0;

    for (int i=0; i < n_species; i++) {
        trace = trace + J(i, i);
    }

    return trace;
}

/*! \brief Python compatible get Jacobian
*   Computes Jacobian from a list of input species, returning a python list containing the Jacobian
*
*   \param input_y A python list of species values used to compute the Jacobian
*   \return Particle Jacobian as python list
*/
boost::python::list Particle::get_jacobian(boost::python::list input_y)
{
    int n_species = state_init.size();
    
    ublas_vec_t y(n_species);
    for (int i=0; i < n_species; i++) {
        y(i) = boost::python::extract<double>(input_y[i]);
    }

    // Init matrix n_species x n_species
    ublas_mat_t J (n_species, n_species);

    // Not sure why this is necessary
    ublas_vec_t dfdt(n_species);

    // Dummy values
    const double t = 0;

    // Fill jacobian matrix
    m.run_jac(y, J, t, dfdt, part_params, model_ref);

    // Unpack ublas jac into python list
    boost::python::list py_J;
    for (int i = 0; i < n_species; i++) {
        for (int j = 0; j < n_species; j++) {
            double val = J(i, j);
            py_J.append(val);
        }
    }

    return py_J;
}

/*! \brief Python compatible get Jacobian from end state
*   Computes Jacobian from a the end state of a simulation, returning a python list containing the Jacobian.
*
*   \return Particle end state Jacobian as python list
*/
boost::python::list Particle::get_end_state_jacobian()
{
    int n_species = state_init.size();
    
    ublas_vec_t y(n_species);
    for (int i=0; i < n_species; i++) {
        y(i) = state_vec.back()[i];
    }

    // Init matrix n_species x n_species
    ublas_mat_t J (n_species, n_species);

    // Not sure why this is necessary
    ublas_vec_t dfdt(n_species);

    // Dummy values
    const double t = 0;

    // Fill jacobian matrix
    m.run_jac(y, J, t, dfdt, part_params, model_ref);

    // Unpack ublas jac into python list
    boost::python::list py_J;
    for (int i = 0; i < n_species; i++) {
        for (int j = 0; j < n_species; j++) {
            double val = J(i, j);
            py_J.append(val);
        }
    }

    return py_J;
}

/*! \brief Python compatible get Jacobian from initial state
* Computes Jacobian using the initial state of the simulation, returning it as a python list
* \return Particle Jacobian as python list
*/
boost::python::list Particle::get_init_state_jacobian()
{
    int n_species = state_init.size();
    
    ublas_vec_t y(n_species);
    for (int i=0; i < n_species; i++) {
        y(i) = state_init[i];
    }

    // Init matrix n_species x n_species
    ublas_mat_t J (n_species, n_species);

    // Not sure why this is necessary
    ublas_vec_t dfdt(n_species);

    // Dummy values
    const double t = 0;

    // Fill jacobian matrix
    m.run_jac(y, J, t, dfdt, part_params, model_ref);

    // Unpack ublas jac into python list
    boost::python::list py_J;
    for (int i = 0; i < n_species; i++) {
        for (int j = 0; j < n_species; j++) {
            double val = J(i, j);
            py_J.append(val);
        }
    }
    return py_J;
}

/*! \brief Gets particle simulation states as a Python list 
* Iterates the state vector, returning a flattened boost::python::list of the results.
* The returned python list reshaping, (n_timepoints, n_species).
* \return sol flat python list of states for a simulated particle
*/
boost::python::list Particle::get_state_pylist() {
    boost::python::list sol;
    {
        ScopedGILRelease noGil = ScopedGILRelease();
        for (auto sim_iter = state_vec.begin(); sim_iter != state_vec.end(); ++sim_iter) {
            auto sim_vec = *sim_iter;
           for(int i=0; i < sim_vec.size(); i++) {
                double val = sim_vec[i];
                sol.append(val);
           }
        }

    }
    return sol;
}

/*! \brief Gets vector of states from particle simulation. 
* Returns vector containing the state list of a simulated particle
* \return std::vector<state_type>&state_vec
*/
std::vector<state_type>& Particle::get_state_vec() {
    return state_vec;
}

/*! \brief Set particle distances
* Sets the particle simulation distances
* \param sim_dist a vector containing simulation distances
*/
void Particle::set_distance_vector(std::vector<std::vector<double>> sim_dist) {
    this->sim_distances = sim_dist;
}

/*! \brief Get particle final species values
* Gets the end state of the simulation, returning a python list of species values
* \return end_state python list of species values for the end of the simulation
*/
boost::python::list Particle::get_final_species_values()
{
    int n_species = state_init.size();

    boost::python::list end_state;

    ublas_vec_t y(n_species);
    for (int i=0; i < n_species; i++) {
        end_state.append(state_vec.back()[i]);
    }

    return end_state;
}

/*! \brief Get sum of simulation standard deviations
* Gets the sum of standard deviations for all species in a particle from a given timepoint to
* the end of the simulation.
* \param from_time_point an integer index of where to start calculating the standard deviations from.
* \return end_state python list of species values for the end of the simulation
*/
double Particle::get_sum_stdev(int from_time_point)
{
    DistanceFunctions dist = DistanceFunctions();

    int n_species = state_init.size();
    return dist.get_sum_stdev(state_vec, n_species, from_time_point);
}

/*! \brief Get final timepoint gradients
* Gets the gradients of all species at the final timepoint.
* \return dist.get_all_species_grads(state_vec, n_species) python list containing the gradients of all species at the
* final timepoint
*/
boost::python::list Particle::get_all_grads()
{
    DistanceFunctions dist = DistanceFunctions();
    int n_species = state_init.size();

    return dist.get_all_species_grads(state_vec, n_species);
}