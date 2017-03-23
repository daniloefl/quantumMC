/**
 * \class PathIntegralMC
 *
 * \ingroup quantumMC
 *
 * \brief Implements the path integral MCMC method for a specific potential.
 *        this effectively solves the Schroedinger equation for that Hamiltonian.
 *
 */

#ifndef PATHINTEGRALMC_H
#define PYTHINTEGRALMC_H

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

class PathIntegralMC {
  public:

    /// \brief Constructor with maximum time and number of MCMC steps required as inputs.
    /// \param potential Python function called with a single double parameter that returns the value of the potential at a certain point
    /// \param xmin Minimum value of the X axis.
    /// \param xmax Maximum value of the X axis.
    /// \param dx Delta X.
    /// \param dx Delta t.
    /// \param Nt Time steps.
    /// \param delta Metropolis X step.
    /// \param reqSteps Number of steps required from the MCMC
    PathIntegralMC(boost::python::object potential, boost::python::object derivPotential,
                   double xmin=-5, double xmax=5, double dx=0.1,
                   double dt = 0.1, int Nt = 100, double delta = 1.0, int reqSteps = 4000);

    /// \brief Set minimum x
    /// \param xmin Minimum value of x
    void setXmin(double xmin);

    /// \brief Set maximum x
    /// \param xmax Maximum value of x
    void setXmax(double xmax);

    /// \brief Set delta x
    /// \param dx Delta x
    void setDeltaX(double dx);

    /// \brief Set MCMC delta
    /// \param dx Delta x
    void setDelta(double dx);

    /// \brief Set delta t and number of time steps
    /// \param dt Delta t
    /// \param Nt Number of time steps
    void setDeltaT(double dt, int Nt);

    /// \brief Set number of steps
    /// \param reqSteps number of required MCMC steps
    void setNSteps(int reqSteps);

    /// \brief Returns the value of the potential energy in position r
    /// \param r The position where to calculate the potential energy
    double V(double r);

    /// \brief Returns the value of the potential energy derivative in position r
    /// \param r The position where to calculate the potential energy
    double dVdx(double r);

    /// \brief Thermalise and run real MCMC
    void run();

    /// \brief Calculate the ground state energy value
    /// \return Ground state energy value
    double eMean();

    /// \brief Calculate the ground state energy error
    /// \return Ground state energy error
    double eError();

    /// \brief Returns the wave function as a Python list
    /// \return List of tuples (x, psi(x))
    boost::python::list getPsi();

    boost::python::list getEnergy();

  private:
    /// Potential function in Python
    boost::python::object m_potential;

    /// Potential derivative in Python
    boost::python::object m_derivPotential;

    /// Solution of the diffusion equation: the wave function
    std::vector<double> m_psi;

    /// Position of the walkers
    std::vector<double>  m_x;

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

    int m_Nt; /// Number of time steps

    /// Number of required MCMC steps
    int m_reqSteps;

    double m_delta; /// step for Metropolis

    /// \brief Calculate the normalisation of the psi variable
    /// \return Normalisation of psi
    double psiNorm();


    /// \brief Implement an MCMC step.
    /// \param x_n new x value found
    bool step(double &x_n);

    /// \brief Run MCMC for the predetermined set of steps
    void MC();

    /// \brief Run MCMC for 20% of the required steps to thermalise
    void thermalise();

    /// \brief Clean up all variables to get ready for MCMC simulation
    void clean();

    /// Helper variables to generate random numbers
    std::mt19937_64 rEngine;
};

#endif

