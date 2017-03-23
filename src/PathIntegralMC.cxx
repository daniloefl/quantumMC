#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <tuple>
#include <cstdlib>
#include <boost/range/irange.hpp>
#include <boost/python/exec.hpp>
#include <boost/python/extract.hpp>
#include <sstream>
#include <string>
#include <iomanip>

#include <Python.h>

#include "PathIntegralMC.h"

using namespace boost;
using namespace std;

PathIntegralMC::PathIntegralMC(boost::python::object potential, boost::python::object derivPotential, 
                               double xmin, double xmax, double dx,
                               double dt, int Nt, double delta, int reqSteps)
  : rEngine(std::random_device()()) {
  Py_Initialize();

  // Python function holding the potential
  m_potential = potential;
  m_derivPotential = derivPotential;

  // set number of required steps
  m_reqSteps = reqSteps;

  // set axis range
  m_xmin = xmin;
  m_xmax = xmax;
  m_dx = dx;

  // Metropolis step size
  m_delta = delta;

  // set time step
  m_dt = dt;

  // number of time steps
  m_Nt = Nt;

  // set initial energy values
  m_sumE = 0;
  m_sumE2 = 0;

  // set initial size of the wave function
  m_psi.resize(int((m_xmax - m_xmin)/m_dx));
  // position and state of the walkers
  m_x.resize(m_Nt);

  // now create the x position of the walkers
  uniform_real_distribution<double> rFlat5(m_xmin, m_xmax);
  for (int i : irange<int>(0, m_Nt)) {
    // set their positions randomly
    m_x[i] = rFlat5(rEngine);
  }

  // clean it up
  clean();
}


void PathIntegralMC::setXmin(double xmin) {
  m_xmin = xmin;
}

void PathIntegralMC::setXmax(double xmax) {
  m_xmax = xmax;
}

void PathIntegralMC::setDeltaX(double dx) {
  m_dx = dx;
}

void PathIntegralMC::setDelta(double dx) {
  m_delta = dx;
}

void PathIntegralMC::setDeltaT(double dt, int Nt) {
  m_dt = dt;
  m_Nt = Nt;
}


void PathIntegralMC::setNSteps(int reqSteps) {
  m_reqSteps = reqSteps;
}

void PathIntegralMC::clean() {
  m_nMCSteps = 0; // number of taken MC steps
  // reset energy estimate
  m_sumE = 0;
  m_sumE2 = 0;
  // reset wave function as it will be "histogrammed"
  for (auto &x : m_psi) { x = 0.0; }
}

double PathIntegralMC::V(double r) {
  return boost::python::extract<double>(m_potential(r));
}

double PathIntegralMC::dVdx(double r) {
  return boost::python::extract<double>(m_derivPotential(r));
}

bool PathIntegralMC::step(double &x_n) {
  uniform_int_distribution<int> rFlatInt(0, m_Nt-1);
  // choose random time slice:
  int t = rFlatInt(rEngine);

  int t_prev = t-1;
  if (t_prev < 0) t_prev = m_Nt-1;
  int t_next = t+1;
  if (t_next > m_Nt-1) t_next = 0;

  uniform_real_distribution<double> rFlat(-1, 1);
  // shift walker in time t
  // P = exp( - m (x_n - x_{n-1})^2/dt ) incorporates the kinetic energy
  // and it is simulated by the random displacement below
  x_n = m_x[t] + m_delta*rFlat(rEngine); // shifts walker uniformly
  // the equation above is the conditional probability of moving from one
  // step to the next: G(x'|x) = G(x|x') = x + delta*flat(-1,1)

  // now check if we should accept it
  // to generate paths with probability distribution given by
  // the density function, we need to use it in the balance probability equation
  // trial is accepted with probability A
  // A(x'|x) = min (1, P(x') g(x|x') / [P(x) g(x'|x)] )
  double deltaE = 0;  // change in the energy
  deltaE += V(x_n) - V(m_x[t]); // change in potential energy
  deltaE += 0.5*std::pow((m_x[t_next] - x_n)/m_dt, 2) - 0.5*std::pow((m_x[t_next] - m_x[t])/m_dt, 2); // change in kinetic energy
  deltaE += 0.5*std::pow((x_n - m_x[t_prev])/m_dt, 2) - 0.5*std::pow((m_x[t] - m_x[t_prev])/m_dt, 2); // change in kinetic energy
  double A = std::exp(-m_dt*deltaE);
  if (A > 1) A = 1;
  uniform_real_distribution<double> rFlatA(0, 1);
  if (A > rFlatA(rEngine)) { // accept
    m_x[t] = x_n;
    return true;
  }
  x_n = m_x[t];
  return false;
}

void PathIntegralMC::MC() {
  for (int i : irange<int>(0, m_Nt)) {
    double x_n = 0;
    bool accepted = step(x_n); // make Nt MCMC steps

    double E0 = V(x_n) + 0.5*x_n*dVdx(x_n); // kinetic energy from the Virial theorem
    // calculate average energy and average energy^2
    m_sumE += E0;
    m_sumE2 += E0*E0;

    // find the index of the walker in psi
    int k = int((x_n-m_xmin)/(m_xmax-m_xmin)*m_psi.size());
    // increment psi in that index to histogram it
    if (k < m_psi.size() && k >= 0)
      m_psi[k] += 1;
  }
  m_nMCSteps++; // end of this MCMC step
}

void PathIntegralMC::thermalise() {
  // just do 20% of the requested steps to initialise walkers to something
  for (int i : irange<int>(0, int(0.2*m_reqSteps))) {
    MC();
  }
}

void PathIntegralMC::run() {
  clean(); // clean up results
  thermalise(); // initialise walkers with something close to the final distribution
  clean(); // clean up results but keep walkers in the thermalised positions
  for (int i : irange<int>(0, m_reqSteps)) { // now do the requested MCMC steps
    MC();
  }
}

double PathIntegralMC::eMean() {
  // get mean energy
  return m_sumE/double(m_nMCSteps*m_Nt);
}

double PathIntegralMC::eError() {
  // get sqrt(variance) of the energy estimate
  return sqrt(m_sumE2/double(m_nMCSteps*m_Nt) - pow(m_sumE/double(m_nMCSteps*m_Nt), 2))/sqrt(m_nMCSteps*m_Nt);
}

double PathIntegralMC::psiNorm() {
  double n = 0;
  // get the sum of |psi|^2 delta x to estimate integral |psi|^2 dx
  for (int i : irange<int>(0, (int) m_psi.size())) n += pow(m_psi[i], 2)*m_dx;
  return n;
}

python::list PathIntegralMC::getPsi() {
  double norm = psiNorm();
  python::list l;
  python::list x;
  python::list psi;
  for (int i : irange<int>(0, (int) m_psi.size())) {
    x.append(i*m_dx + m_xmin);
    psi.append(m_psi[i]/sqrt(norm));
  }
  l.append(x);
  l.append(psi);
  return l;
}

python::list PathIntegralMC::getEnergy() {
  python::list l;
  l.append(eMean());
  l.append(eError());
  return l;
}

