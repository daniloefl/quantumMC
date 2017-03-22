#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>

using namespace std;

random_device r;
mt19937_64 rEngine(r());
normal_distribution<double> rGaus(0.0, 1.0);


// trial function here
class WF {
  public:
    //double m_alpha; // to be chosen
    vector<double> m_alpha; // to be chosen

    vector<double> m_psi2;  // histogram to be filled
    vector<double> m_x;     // walkers

    double m_sumE;          // sum of energy from accepted MC
    double m_sumE2;
    double m_nAccept;       // number of accepted iteractions
    double m_nTotal;        // number of total iteractions
    double m_nMCSteps;      // number of total MC steps

    // inputs
    double m_xmin;
    double m_xmax;
    double m_dx;
    int m_N;
    double m_sigma;
    int m_thermSteps;

    void clean() {
      m_nTotal = 0;
      m_nAccept = 0;
      m_nMCSteps = 0;
      m_sumE = 0;
      m_sumE2 = 0;
      for (auto &x : m_psi2) { x = 0; }
    }

    WF(int N = 300) {
      // could be inputs
      m_N = N;
      m_xmin = -10;
      m_xmax = +10;
      m_dx = 0.1;
      m_sigma = 1.0;
      m_thermSteps = 1000;

      m_sumE = 0;
      m_sumE2 = 0;

      m_psi2.resize(int((m_xmax - m_xmin)/m_dx));
      m_x.resize(m_N);
      uniform_real_distribution<double> rFlat5(-0.5, 0.5);
      for (auto &x : m_x) { x = rFlat5(rEngine); }
      m_alpha.resize(1);

      clean();
    }

    double psiT(double x) {
      return exp(-m_alpha[0]*x*x); // harmonic oscillator
      //double r = 0;
      //for (int n = 0; n < m_alpha.size(); ++n) {
      //  r += m_alpha[n]*pow(x, n);
      //}
      //return r;
    }

    // H |psi> = E_L |psi>
    // H = -1/2 d^2/dx^2 + 1/2 x^2  (m = 1, omega = 1)
    // E_L|psi> = alpha exp(-alpha x^2) - 2*alpha^2 x^2 exp(-alpha x^2) + 0.5 * x^2 exp(-alpha x^2)
    // E_L = alpha - 2 alpha^2 x^2 + 0.5 x^2
    double eL(double x) {
      return m_alpha[0] - 2.0*pow(m_alpha[0]*x, 2) + 0.5*pow(x, 2);

      //double r = 0;
      //for (int n = 2; n < m_alpha.size(); ++n) {
      //  r += -0.5*n*(n-1)*m_alpha[n]*pow(x, n-2);
      //}
      //return r + 0.5*pow(x, 2);
    }

    // gives probability distribution of |psi(x)|^2
    // used to generate numbers according to the trial distribution
    double rho(double x) {
      return pow(psiT(x), 2);
    }

    void step() {
      // choose a random element
      uniform_int_distribution<> rIntN(0, m_N-1);
      uniform_real_distribution<double> rFlat(0, 1.0);
      int n = rIntN(rEngine);
      // shift walker
      double xTrial = m_x[n] + m_sigma*rGaus(rEngine);
      // MC step:
      // if p(new x)/p(old x) > flat number between 0 and 1,
      // then accept it, otherwise reject it
      // this generates x[n] with distribution given by p
      if (rho(xTrial)/rho(m_x[n]) > rFlat(rEngine)) {
        m_x[n] = xTrial;
        m_nAccept++;
      }
      m_nTotal++;

      m_sumE += eL(m_x[n]);
      m_sumE2 += pow(eL(m_x[n]), 2);

      int i = (int) ((m_x[n] - m_xmin)/m_dx);
      if (i < m_psi2.size()) m_psi2[i] += 1;
    }

    void MC() {
      for (int i = 0; i < m_N; ++i) step();
      m_nMCSteps++;
    }

    void thermalise() {
      int adjSteps = int(0.1*m_thermSteps);
      for (int i = 0; i < m_thermSteps; ++i) {
        MC();
        if (i % adjSteps == 0) {
          m_sigma *= m_nAccept/(0.5*m_N*adjSteps);
          m_nAccept = 0;
        }
      }
    }

    double eMean() {
      return m_sumE/double(m_N)/double(m_nMCSteps);
    }

    double eError() {
      return sqrt(m_sumE2/double(m_N)/double(m_nMCSteps) - pow(m_sumE/double(m_N)/double(m_nMCSteps), 2))/sqrt(m_N*m_nMCSteps);
    }

    double psiNorm() {
      double n = 0;
      for (auto &i : m_psi2) n += i*m_dx;
      return n;
    }

    vector<double> &alpha() {
      return m_alpha;
    }


    void write(const string &f) {
      ofstream ff(f.c_str());
      //ff << "# $|\\psi|^2 = ";
      //for (int n = 0; n < m_alpha.size(); ++n) {
      //  ff << m_alpha[n] << " x^" << n;
      //  if (n < m_alpha.size() - 1)
      //    ff << " + ";
      //}
      //ff << "$" << endl;
      ff << "# $|\\psi|^2 = exp(-"<< m_alpha[0] << "x^2)$" << endl;
      ff << "# $E = " << eMean() << " \\pm " << eError() << "$" << endl;
      ff << "# x          psi^2" << endl;
      double norm = psiNorm();
      for (int i = 0; i < m_psi2.size(); ++i) {
        ff << i*m_dx + m_xmin << "    " << m_psi2[i]/norm << endl;
      }
      ff.close();
    }
};

int main(int argc, char **argv) {

  WF psi;
  psi.alpha().resize(1);
  
  psi.thermalise();
  psi.clean();


  double dalpha = 0.05;
  double minalpha = 0.3;
  int Nalpha = 12;

  double minE = 99;
  double bestAlpha = 10;
  ofstream ff("testE.dat");
  ff << "# alpha    E    dE" << endl;
  for (int i = 0; i < Nalpha+1; ++i) {
    psi.alpha()[0] = minalpha + i*dalpha;

    psi.clean();
    int MCSteps = 10000;
    for (int i = 0; i < MCSteps; ++i)
      psi.MC();
    if (psi.eMean() < minE) {
      minE = psi.eMean();
      bestAlpha = psi.alpha()[0];
      psi.write("psi.dat");
    }
    ff << psi.alpha()[0] << "    "<< psi.eMean() << "    " << psi.eError() << endl;
  }
  ff.close();

  return 0;
}

