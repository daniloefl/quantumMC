/**
 * \class QuantumMC
 *
 * \ingroup quantumMC
 *
 * \brief Implements diffusion MCMC method for a specific potential in 3D.
 *
 */

#ifndef QUANTUMMC_H
#define QUANTUMMC_H

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

#define DIM 1
typedef std::array<double, DIM> pos;

class QuantumMC {
  public:

    /// \brief Constructor with maximum time and number of MCMC steps required as inputs.
    /// \param NT
    /// \param reqSteps Number of steps required from the MCMC
    WF(int NT = 10000, int reqSteps = 4000);

    /// \brief Returns the value of the potential energy in position r
    /// \param r The position where to calculate the potential energy
    double V(pos r);

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

  private:
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

};

#endif

