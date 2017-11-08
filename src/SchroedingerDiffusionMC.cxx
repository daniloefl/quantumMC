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

#include "SchroedingerDiffusionMC.h"

using namespace boost;
using namespace std;

void SchroedingerDiffusionMC::setImportanceSampling(bool active, boost::python::object guidingWF, boost::python::object localEnergy, boost::python::object quantumForce) {
  m_importanceSampling = active;

  // Python function holding the guiding wave function
  m_guidingWF = guidingWF;

  // Python function holding the local energy
  m_localEnergy = localEnergy;

  // Python function holding the quantum force
  m_quantumForce = quantumForce;

}

SchroedingerDiffusionMC::SchroedingerDiffusionMC(boost::python::object potential,
                                                 double xmin, double xmax, double dx,
                                                 int NT, int reqSteps)
  : rEngine(std::random_device()()) {

  Py_Initialize();

  m_logGrid = false;

  // Python function holding the potential
  m_potential = potential;

  // set number of required steps
  m_reqSteps = reqSteps;

  m_EL = 0;
  m_countEL = 0;
  m_sumEL = 0;
  m_sumEL2 = 0;

  // number of surviving walkers
  m_NT = NT;

  // number of walkers
  m_N = m_NT;

  // set axes range
  m_xmin = xmin;
  m_xmax = xmax;
  m_dx = dx;

  // set time step
  m_dt = 0.1;

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
  uniform_real_distribution<double> rFlat5(m_xmin, m_xmax);
  for (int i : irange<int>(0, m_N)) {
    m_x[i] = pos();
    // set their positions randomly
    m_x[i][0] = rFlat5(rEngine);
    m_alive[i] = true; // all are alive
  }

  // clean it up
  clean();
}

void SchroedingerDiffusionMC::setTimeStep(double dt) {
  m_dt = dt;
}

void SchroedingerDiffusionMC::setXmin(double xmin) {
  m_xmin = xmin;
}

void SchroedingerDiffusionMC::setXmax(double xmax) {
  m_xmax = xmax;
}

void SchroedingerDiffusionMC::setDeltaX(double dx) {
  m_dx = dx;
}

void SchroedingerDiffusionMC::logGrid(bool logGrid) {
  m_logGrid = logGrid;
}


void SchroedingerDiffusionMC::setNSteps(int reqSteps) {
  m_reqSteps = reqSteps;
}

void SchroedingerDiffusionMC::setN(int NT) {
  m_NT = NT;
  m_N = NT;
}

void SchroedingerDiffusionMC::clean() {
  m_nMCSteps = 0; // number of taken MC steps
  // reset energy estimate
  m_sumE = 0;
  m_sumE2 = 0;
  m_EL = 0;
  m_countEL = 0;
  m_sumEL = 0;
  m_sumEL2 = 0;
  // reset wave function as it will be "histogrammed"
  for (auto &x : m_psi) { x = 0.0; }
}

double SchroedingerDiffusionMC::V(pos r) {
  //return 0.5*(pow(r[0], 2)); // + pow(r[1], 2) + pow(r[2], 2));

  //std::stringstream ss;
  //ss << m_potentialName << "(" << setprecision(17) << r[0] << ")";
  //boost::python::object V_obj = boost::python::eval(ss.str().c_str());
  //double V_double = boost::python::extract<double>(V_obj);

  return boost::python::extract<double>(m_potential(r[0]));
}

double SchroedingerDiffusionMC::V(double r) {
  return boost::python::extract<double>(m_potential(r));
}

bool SchroedingerDiffusionMC::step(int n) {
  uniform_real_distribution<double> rFlat(0, 1);
  // shift walker
  // psi = prod dx_j [ prod W P(x_n, x_{n-1}) ] psi(0)
  // Weight function W = exp(-dt*(V - E))

  // P = exp( - m (x_n - x_{n-1})^2/dt ) incorporates the kinetic energy
  // and it is simulated by the random displacement below
  normal_distribution<double> rGaus(0.0, 1.0);
  double xo = m_x[n][0];
  double xn = xo + sqrt(m_dt)*rGaus(rEngine);
  m_x[n][0] = xn;

  double ELo = 0;
  double EL = 0;
  if (m_logGrid) {
    EL = V(std::exp(xn));
    ELo = V(std::exp(xo));
  } else {
    EL = V(xn);
    ELo = V(xo);
  }

  double Wb = 1;
  double We = 1;
  if (m_logGrid) {
    Wb *= std::exp(-m_dt*0.5*(std::pow(std::exp(xn), 2) + std::pow(std::exp(xo), 2))*(0.5*(EL+ELo) - m_E0));
    We *= std::exp(-m_dt*(0.5*(EL+ELo) - m_E0));
  } else {
    Wb *= std::exp(-m_dt*(0.5*(EL+ELo) - m_E0));
    We = Wb;
  }
  m_EL += EL*We;
  m_countEL += We;
  int survivors = int(Wb);
  if (Wb - survivors > rFlat(rEngine)) {
    survivors++;
  }
  // survivors = W with probability frac(W) = W - int(W)
  // and survivors = W+1 with probability 1 - frac(W)
  // now we have an integer number of new particles
  // make a new walker
  if (survivors-1 > 0) {
    m_x.resize(m_N+(survivors-1));
    m_alive.resize(m_N+(survivors-1));
    for (int i : irange<int>(0, survivors - 1)) {
      m_x[m_N][0] = m_x[n][0];
  
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
  return true;
}

bool SchroedingerDiffusionMC::stepImportanceSampling(int n) {
  uniform_real_distribution<double> rFlat(0, 1);
  // shift walker
  // psi = prod dx_j [ prod W P(x_n, x_{n-1}) ] psi(0)
  // Weight function W = exp(-dt*(V - E))

  // P = exp( - m (x_n - x_{n-1})^2/dt ) incorporates the kinetic energy
  // and it is simulated by the random displacement below
  normal_distribution<double> rGaus(0.0, 1.0);
  double xo = m_x[n][0];

  double Fo = 0;
  if (m_logGrid) Fo = quantumForce(std::exp(xo));
  else Fo = quantumForce(xo);

  double xn = xo + sqrt(m_dt)*rGaus(rEngine) + m_dt*Fo;
  double Fn = 0;
  if (m_logGrid) Fn = quantumForce(std::exp(xn));
  else Fn = quantumForce(xn);


  double ELo = 0;
  double EL = 0;
  if (m_logGrid) {
    EL = localEnergy(std::exp(xn));
    ELo = localEnergy(std::exp(xo));
  } else {
    EL = localEnergy(xn);
    ELo = localEnergy(xo);
  }

  double W = 1;
  if (m_logGrid) {
    W *= std::exp(-std::pow(std::exp(xo) - std::exp(xn) - m_dt*Fn, 2)/(2*m_dt));
    W /= std::exp(-std::pow(std::exp(xn) - std::exp(xo) - m_dt*Fo, 2)/(2*m_dt));
    W *= std::pow(guidingWF(std::exp(xn)), 2);
    W /= std::pow(guidingWF(std::exp(xo)), 2);
    W *= std::exp(-m_dt*(std::pow(std::exp(xn), 2))*(0.5*(EL+ELo) - m_E0));
    W /= std::exp(-m_dt*(std::pow(std::exp(xo), 2))*(0.5*(EL+ELo) - m_E0));
  } else {
    W *= std::exp(-std::pow(xo - xn - m_dt*Fn, 2)/(2*m_dt));
    W /= std::exp(-std::pow(xn - xo - m_dt*Fo, 2)/(2*m_dt));
    W *= std::pow(guidingWF(xn), 2);
    W /= std::pow(guidingWF(xo), 2);
  }

  double A = 1;
  if (W < 1) A = W;

  if (rFlat(rEngine) < A) { // probability A
    // accept

    m_x[n][0] = xn;

    double Wb = 1;
    double We = 1;
    if (m_logGrid) {
      //Wb = std::exp(-m_dt*(std::pow(std::exp(xn), 2))*(0.5*(EL+ELo) - m_E0));
      Wb = std::exp(-m_dt*0.5*(std::pow(std::exp(xn), 2) + std::pow(std::exp(xo), 2))*(0.5*(EL+ELo) - m_E0));
      We = Wb;
    } else {
      Wb = std::exp(-m_dt*(0.5*(EL+ELo) - m_E0));
      We = Wb;
    }
    m_EL += EL*We;
    m_countEL += We;
    int survivors = int(Wb);
    if (Wb - survivors > rFlat(rEngine)) {
      survivors++;
    }
    // survivors = W with probability frac(W) = W - int(W)
    // and survivors = W+1 with probability 1 - frac(W)
    // now we have an integer number of new particles
    // make a new walker
    if (survivors-1 > 0) {
      m_x.resize(m_N+(survivors-1));
      m_alive.resize(m_N+(survivors-1));
      for (int i : irange<int>(0, survivors - 1)) {
        m_x[m_N][0] = m_x[n][0];
  
        m_alive[m_N] = true;
        m_N++;
      }
    }
    // kill the walker
    if (survivors == 0)
      m_alive[n] = false;
    // too many walkers ...
    if (m_N > 100000) {
      cout << "Too many nodes: " << m_N << std::endl;
      exit(-1);
    }
    return true;
  }
  return false;
}

void SchroedingerDiffusionMC::MC() {
  int N0 = m_N; // original number of walkers
  m_EL = 0;
  double accep = 0;
  if (m_importanceSampling) {
    for (int i : irange<int>(0, N0)) {
      if (stepImportanceSampling(i)) {
        accep += 1; // make N0 MCMC steps
      }
    }
  } else {
    for (int i : irange<int>(0, N0)) {
      if (step(i)) {
        accep += 1; // make N0 MCMC steps
      }
    }
  }

  // now count the number of alive walkers
  // this also cleans up the vectors removing the dead walkers from the list
  int j = 0;
  for (int i : irange<int>(0, m_N)) {
    if (m_alive[i] && i != j) {
      m_x[j] = m_x[i];
      m_alive[j] = m_alive[i];
    }
    if (m_alive[i]) ++j;
  }
  m_N = j; // reset number of alive walkers
  // https://www.thphys.uni-heidelberg.de/~wetzel/qmc2006/KOSZ96.pdf
  // page 3
  // the energy must be set to go in the direction of the
  // ground state energy, otherwise the wave function vanishes or diverges
  // see end of page 5, beginning of page 6
  // change in number of walkers implies change in energy
  // setting alpha = 0.1
  m_E0 += log(double(m_NT)/double(m_N))*0.1;
  //m_E0 += 1.0/m_dt*(1.0 - ((double) m_N)/((double) N0));
  // calculate average energy and average energy^2
  m_sumE += m_E0;
  m_sumE2 += m_E0*m_E0;
  m_sumEL += m_EL;
  m_sumEL2 += m_EL*m_EL;
  // finally, histogram in this step for each walker
  for (int i : irange<int>(0, m_N)) {
    // find out the position of each walker
    double x = m_x[i][0];
    // find the index of the walker in psi
    int k = int((x-m_xmin)/(m_xmax-m_xmin)*m_psi.size());
    // increment psi in that index to histogram it
    if (k < m_psi.size() && k > 0)
      m_psi[k] += 1;
  }
  m_nMCSteps++; // end of this MCMC step
}

void SchroedingerDiffusionMC::thermalise() {
  // just do 20% of the requested steps to initialise walkers to something
  std::cout << "Thermalisation step using 40\% of the requested time range." << std::endl;
  std::cout << std::setw(40) << "Time step size: " << std::setw(10) << std::setprecision(5) << m_dt << std::endl;
  std::cout << std::setw(40) << "Thermalisation total time goal: " << std::setw(10) << std::setprecision(5) << 0.4*m_reqSteps*m_dt << std::endl;
  std::cout << std::setw(10) << "Time" << " " << std::setw(25) << "Average local energy" << " " << std::setw(25) << "Average trial energy" << std::endl;
  for (int i : irange<int>(0, int(0.4*m_reqSteps))) {
    if (i > 0 && i % 1000 == 0) {
      std::cout << std::setw(10) << std::setprecision(5) << i*m_dt << " " << std::setw(10) << std::setprecision(8) << eLMean() << " +/- " << std::setw(10) << std::setprecision(8) << eLError() << " " << std::setw(10) << std::setprecision(8) << eMean() << " +/- " << std::setw(10) << std::setprecision(8) << eError() << std::endl;
    }
    MC();
  }
}

void SchroedingerDiffusionMC::run() {
  // set initial size of the wave function
  m_psi.resize(int((m_xmax - m_xmin)/m_dx));
  // position and state of the walkers
  m_x.resize(m_N);
  m_alive.resize(m_N);

  // now create the x position of the walkers
  uniform_real_distribution<double> rFlat5(m_xmin, m_xmax);
  for (int i : irange<int>(0, m_N)) {
    m_x[i] = pos();
    // set their positions randomly
    m_x[i][0] = rFlat5(rEngine);
    m_alive[i] = true; // all are alive
  }

  // clean it up
  clean(); // clean up results
  thermalise(); // initialise walkers with something close to the final distribution
  clean(); // clean up results but keep walkers in the thermalised positions

  std::cout << "Beginning calculation." << std::endl;
  std::cout << std::setw(40) << "Time step size: " << std::setw(10) << std::setprecision(5) << m_dt << std::endl;
  std::cout << std::setw(40) << "Requested total time: " << std::setw(10) << std::setprecision(5) << m_reqSteps*m_dt << std::endl;

  std::cout << std::setw(10) << "Time" << " " << std::setw(25) << "Average local energy" << " " << std::setw(25) << "Average trial energy" << std::endl;
  for (int i : irange<int>(0, m_reqSteps)) { // now do the requested MCMC steps
    if (i > 0 && i % 1000 == 0) {
      std::cout << std::setw(10) << std::setprecision(5) << i*m_dt << " " << std::setw(10) << std::setprecision(8) << eLMean() << " +/- " << std::setw(10) << std::setprecision(8) << eLError() << " " << std::setw(10) << std::setprecision(8) << eMean() << " +/- " << std::setw(10) << std::setprecision(8) << eError() << std::endl;
    }
    MC();
  }
}

double SchroedingerDiffusionMC::eMean() {
  // get mean energy
  return m_sumE/double(m_nMCSteps);
}

double SchroedingerDiffusionMC::eError() {
  // get sqrt(variance) of the energy estimate
  return sqrt((m_sumE2/double(m_nMCSteps) - pow(m_sumE/double(m_nMCSteps), 2))/((double) (m_nMCSteps - 1)));
}

double SchroedingerDiffusionMC::eLMean() {
  // get mean energy
  return m_sumEL/double(m_countEL);
}

double SchroedingerDiffusionMC::eLError() {
  // get sqrt(variance) of the energy estimate
  return sqrt((m_sumEL2/double(m_countEL) - pow(m_sumEL/double(m_countEL), 2))/((double) (m_countEL - 1)));
}

double SchroedingerDiffusionMC::localEnergy(double x) {
  return boost::python::extract<double>(m_localEnergy(x));
}

double SchroedingerDiffusionMC::quantumForce(double x) {
  return boost::python::extract<double>(m_quantumForce(x));
}

double SchroedingerDiffusionMC::guidingWF(double x) {
  return boost::python::extract<double>(m_guidingWF(x));
}

double SchroedingerDiffusionMC::psiNorm() {
  double n = 0;
  // get the sum of |psi|^2 delta x to estimate integral |psi|^2 dx
  if (m_logGrid) {
    for (int i : irange<int>(0, (int) m_psi.size())) n += pow(m_psi[i], 2)*std::pow(std::exp(i*m_dx + m_xmin), 2)*m_dx;
  } else {
    for (int i : irange<int>(0, (int) m_psi.size())) n += pow(m_psi[i], 2)*m_dx;
  }
  return n;
}

void SchroedingerDiffusionMC::write(const string &f) {
  // write out the result
  ofstream ff(f.c_str());
  ff << "# $E = " << eMean() << " \\pm " << eError() << "$" << endl;
  ff << "# r          psi" << endl;
  double norm = psiNorm();
  for (int i : irange<int>(0, (int) m_psi.size())) {
    //ff << i*m_dx + m_xmin << "    " << pow(i*m_dx, DIM-1)*pow(m_psi[i], 2)/norm << endl;
    if (m_logGrid) {
      ff << std::exp(i*m_dx + m_xmin) << "    " << m_psi[i]*std::sqrt(std::exp(i*m_dx + m_xmin))/sqrt(norm) << endl;
    } else {
      ff << i*m_dx + m_xmin << "    " << m_psi[i]/sqrt(norm) << endl;
    }
  }
  ff.close();
}

python::list SchroedingerDiffusionMC::getPsi() {
  double norm = psiNorm();
  python::list l;
  python::list x;
  python::list psi;
  for (int i : irange<int>(0, (int) m_psi.size())) {
    //l.append(python::make_tuple(i*m_dx + m_xmin, m_psi[i]/sqrt(norm)));
    if (m_logGrid) {
      x.append(std::exp(i*m_dx + m_xmin));
      psi.append(m_psi[i]*std::pow(std::exp(i*m_dx + m_xmin), 0.5)/sqrt(norm));
    } else {
      x.append(i*m_dx + m_xmin);
      psi.append(m_psi[i]/sqrt(norm));
    }
  }
  l.append(x);
  l.append(psi);
  return l;
}

python::list SchroedingerDiffusionMC::getEnergy() {
  python::list l;
  l.append(eMean());
  l.append(eError());
  return l;
}

python::list SchroedingerDiffusionMC::getLocalEnergy() {
  python::list l;
  l.append(eLMean());
  l.append(eLError());
  return l;
}

