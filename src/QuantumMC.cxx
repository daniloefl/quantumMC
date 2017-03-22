#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <tuple>
#include <cstdlib>
#include <boost/range/irange.hpp>

#include "QuantumMC.h"

using namespace boost;
using namespace std;

random_device r;
mt19937_64 rEngine(r());
normal_distribution<double> rGaus(0.0, 1.0);

QuantumMC::QuantumMC(int NT, int reqSteps) {
  // set number of required steps
  m_reqSteps = reqSteps;

  // number of surviving walkers
  m_NT = NT;

  // number of walkers
  m_N = m_NT;

  // set axes range
  m_xmin = 0.0;
  m_xmax = 10.0;
  m_dx = 0.1;

  // set time step
  m_dt = 0.05;

  // set initial energy values
  m_E0 = 0;
  m_sumE = 0;
  m_sumE2 = 0;

  // set initial size of the wave function
  m_psi.resize(int((m_xmax - m_xmin)/m_dx));
  // position and state of the walkers
  m_x.resize(m_N);
  m_alive.resize(m_N);

  // now create the x position of the walkers
  uniform_real_distribution<double> rFlat5(-0.5, 0.5);
  for (int i : irange(0, m_N)) {
    m_x[i] = pos();
    // set their positions randomly close to zero
    for (int d : irange(0, DIM)) {
      m_x[i][d] = rFlat5(rEngine);
    }
    m_alive[i] = true; // all are alive
  }

  // clean it up
  clean();
}

void QuantumMC::clean() {
  m_nMCSteps = 0; // number of taken MC steps
  // reset energy estimate
  m_sumE = 0;
  m_sumE2 = 0;
  // reset wave function as it will be "histogrammed"
  for (auto &x : m_psi) { x = 0.0; }
}

double QuantumMC::V(pos r) {
  return 0.5*(pow(r[0], 2)); // + pow(r[1], 2) + pow(r[2], 2));
}

void QuantumMC::step(int n) {
  uniform_real_distribution<double> rFlat(0, 1);
  // shift walker
  for (int i : irange(0, DIM)) {
    m_x[n][i] += sqrt(m_dt)*rGaus(rEngine); // shifts walker by a gaussian with width sqrt(delta t): diffusion!
  }
  // change in effective potential energy V - E causes the walker to die with certain probability
  double q = exp(-m_dt*(V(m_x[n]) - m_E0));
  int survivors = int(q);
  if (q - survivors > rFlat(rEngine)) {
    survivors++;
  }
  // make a new walker
  if (survivors-1 > 0) {
    m_x.resize(m_N+(survivors-1));
    m_alive.resize(m_N+(survivors-1));
    for (int i : irange(0, survivors - 1)) {
      for (int d : irange(0, DIM))
        m_x[m_N][d] = m_x[n][d];

      m_alive[m_N] = true;
      m_N++;
    }
  }
  // kill the walker
  if (survivors == 0)
    m_alive[n] = false;
  // too many walkers ...
  if (m_N > 500000) {
    cout << "Too many nodes: " << m_N << std::endl;
    exit(-1);
  }
}

void QuantumMC::MC() {
  int N0 = m_N; // original number of walkers
  for (int i : irange(0, N0)) step(i); // make N0 MCMC steps

  // now count the number of alive walkers
  // this also cleans up the vectors removing the dead walkers from the list
  int j = 0;
  for (int i : irange(0, m_N)) {
    if (m_alive[i] && i != j) {
      m_x[j] = m_x[i];
      m_alive[j] = m_alive[i];
    }
    if (m_alive[i]) ++j;
  }
  m_N = j; // reset number of alive walkers
  // change in number of walkers implies change in energy
  // setting alpha = 0.1
  m_E0 += log(double(m_NT)/double(m_N))/10.0;
  // calculate average energy and average energy^2
  m_sumE += m_E0;
  m_sumE2 += m_E0*m_E0;
  // finally, histogram in this step for each walker
  for (int i : irange(0, m_N)) {
    // find out the position of each walker
    double r2 = 0;
    for (int d : irange(0, DIM)) {
      r2 += pow(m_x[i][d], 2);
    }
    // find the index of the walker in psi
    int k = int(sqrt(r2)/m_xmax*m_psi.size());
    // increment psi in that index to histogram it
    if (k < m_psi.size())
      m_psi[k] += 1;
  }
  m_nMCSteps++; // end of this MCMC step
}

void QuantumMC::thermalise() {
  // just do 20% of the requested steps to initialise walkers to something
  for (int i : irange(0, int(0.2*m_reqSteps))) {
    MC();
  }
}

void QuantumMC::run() {
  thermalise(); // initialise walkers with something close to the final distribution
  clean(); // clean up results but keep walkers in the thermalised positions
  for (int i : irange(0, m_reqSteps)) { // now do the requested MCMC steps
    MC();
  }
}

double QuantumMC::eMean() {
  // get mean energy
  return m_sumE/double(m_nMCSteps);
}

double QuantumMC::eError() {
  // get sqrt(variance) of the energy estimate
  return sqrt(m_sumE2/double(m_nMCSteps) - pow(m_sumE/double(m_nMCSteps), 2))/sqrt(m_nMCSteps);
}

double QuantumMC::psiNorm() {
  double n = 0;
  // get the sum of |psi|^2 delta x to estimate integral |psi|^2 dx
  for (int i : irange(0, m_psi.size())) n += pow(i*m_dx, DIM-1)*pow(m_psi[i], 2)*m_dx;
  return n;
}

void QuantumMC::write(const string &f) {
  // write out the result
  ofstream ff(f.c_str());
  ff << "# $E = " << eMean() << " \\pm " << eError() << "$" << endl;
  ff << "# r          psi" << endl;
  double norm = psiNorm();
  for (int i : irange(0, m_psi.size())) {
    //ff << i*m_dx + m_xmin << "    " << pow(i*m_dx, DIM-1)*pow(m_psi[i], 2)/norm << endl;
    ff << i*m_dx + m_xmin << "    " << m_psi[i]/sqrt(norm) << endl;
  }
  ff.close();
}

