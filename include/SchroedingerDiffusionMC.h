/**
 * \class SchroedingerDiffusionMC
 *
 * \ingroup quantumMC
 *
 * \brief Implements diffusion MCMC method for a specific potential in 3D.
 *
 */

#ifndef SCHROEDINGERDIFFUSIONMC_H
#define SCHROEDINGERDIFFUSIONMC_H

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <tuple>
#include <cstdlib>
#include <boost/range/irange.hpp>
#include <string>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

typedef std::array<double, 1> pos;

class SchroedingerDiffusionMC {
  public:

    /// \brief Constructor with maximum time and number of MCMC steps required as inputs.
    /// \param potential Python function called with a single double parameter that returns the value of the potential at a certain point
    /// \param xmin Minimum value of the X axis.
    /// \param xmax Maximum value of the X axis.
    /// \param dx Delta X.
    /// \param NT Number of walkers
    /// \param reqSteps Number of steps required from the MCMC
    SchroedingerDiffusionMC(boost::python::object potential, double xmin=-5, double xmax=5, double dx=0.1, int NT = 300, int reqSteps = 4000);

    /// \brief Set minimum x
    /// \param xmin Minimum value of x
    void setXmin(double xmin);

    /// \brief Set maximum x
    /// \param xmax Maximum value of x
    void setXmax(double xmax);

    /// \brief Set delta x
    /// \param dx Delta x
    void setDeltaX(double dx);

    /// \brief Set logarithmic Grid.
    /// \param logGrid Whether the Grid is to be logarithmic
    void logGrid(bool logGrid = true);

    /// \brief Set number of steps
    /// \param reqSteps number of required MCMC steps
    void setNSteps(int reqSteps);

    /// \brief Set number of walkers
    /// \param NT number of walkers
    void setN(int NT);

    /// \brief Returns the value of the potential energy in position r
    /// \param r The position where to calculate the potential energy
    double V(pos r);

    /// \brief Returns the value of the potential energy in position r
    /// \param r The position where to calculate the potential energy
    double V(double r);

    /// \brief Thermalise and run real MCMC
    void run();

    /// \brief Calculate the ground state energy value
    /// \return Ground state energy value
    double eMean();

    /// \brief Calculate the ground state energy error
    /// \return Ground state energy error
    double eError();

    /// \brief Dump resulting distribution in file named by f
    /// \param f Filename to dump result to
    void write(const std::string &f);

    /// \brief Returns the wave function as a Python list
    /// \return List of tuples (x, psi(x))
    boost::python::list getPsi();

    boost::python::list getEnergy();

  private:
    /// Logarithmic grid
    bool m_logGrid;

    /// Potential function in Python and its object
    boost::python::object m_potential;

    /// Solution of the diffusion equation: the wave function
    std::vector<double> m_psi;

    /// Position of the walkers
    std::vector<pos>  m_x;
    /// Whether the walker is alive or died
    std::vector<bool> m_alive;

    /// Ground state energy
    double m_E0;
    /// Sum of energy from accepted simulations
    double m_sumE;
    /// Sum of energy squared
    double m_sumE2;

    /// Number of MCMC steps
    double m_nMCSteps;

    /// Minimum, maximum and direction step size
    double m_xmin;
    double m_xmax;
    double m_dx;

    /// Time step used for the diffusion evolution
    double m_dt;

    /// Number of walkers alive
    int m_N;

    /// Number of initial walkers
    int m_NT;

    /// Number of required MCMC steps
    int m_reqSteps;

    /// \brief Calculate the normalisation of the psi variable
    /// \return Normalisation of psi
    double psiNorm();


    /// \brief Implement an MCMC step on walker n
    /// \param n Walker number to apply the step on
    void step(int n);

    /// \brief Run MCMC for the predetermined set of steps
    void MC();

    /// \brief Run MCMC for 20% of the required steps to thermalise
    void thermalise();

    /// \brief Clean up all variables to get ready for MCMC simulation
    void clean();

    /// Helper variables for random number generation
    std::mt19937_64 rEngine;
};

#endif

