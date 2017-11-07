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

    /// \brief Activate importance sampling
    /// \param active Whether to activate importance sampling.
    /// \param guidingWF Trial wavefunction.
    /// \param localEnergy Value of the hamiltonian applied in the trial wave function divided by the trial wave function
    /// \param quantumForce Value of the derivative of the trial wave function over the trial wave function
    void setImportanceSampling(bool active, boost::python::object guidingWF, boost::python::object localEnergy, boost::python::object quantumForce);

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

    /// \brief Calculate local energy relative to reference function
    /// \param Position
    double localEnergy(double x);

    /// \brief Calculate quantum force
    /// \param Position
    double quantumForce(double x);

    /// \brief Calculate guiding wavefunction
    /// \param Position
    double guidingWF(double x);

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

    /// \brief Calculate the average local energy
    /// \return Average local energy
    double eLMean();

    /// \brief Calculate the average local energy error
    /// \return Average local energy error
    double eLError();

    /// \brief Dump resulting distribution in file named by f
    /// \param f Filename to dump result to
    void write(const std::string &f);

    /// \brief Returns the wave function as a Python list
    /// \return List of tuples (x, psi(x))
    boost::python::list getPsi();

    boost::python::list getEnergy();

    boost::python::list getLocalEnergy();

  private:
    /// Whether to use importance sampling
    bool m_importanceSampling;

    /// Logarithmic grid
    bool m_logGrid;

    /// Potential function in Python and its object
    boost::python::object m_potential;

    /// Guiding wave function
    boost::python::object m_guidingWF;

    /// Local energy function in Python
    boost::python::object m_localEnergy;

    /// Quantum force in Python
    boost::python::object m_quantumForce;

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

    /// Current local energy relative to reference wave function
    double m_EL;

    /// Accummulation of weight factor for EL normalisation
    double m_countEL;

    /// Sum of local energies
    double m_sumEL;

    /// Sum in squares of local energies
    double m_sumEL2;

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
    /// \return Whether the MCMC event was accepted
    bool step(int n);

    /// \brief Implement an MCMC step on walker n using importance sampling
    /// \param n Walker number to apply the step on
    /// \return Whether the MCMC event was accepted
    bool stepImportanceSampling(int n);

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

