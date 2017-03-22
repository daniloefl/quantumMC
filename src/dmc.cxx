#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <tuple>
#include <cstdlib>

using namespace std;

random_device r;
mt19937_64 rEngine(r());
normal_distribution<double> rGaus(0.0, 1.0);

#define DIM 3
typedef array<double, DIM> pos;

class range
{
private:
    int last;
    int iter;

public:
    range(int end):
        last(end),
        iter(0)
    {}

    // Iterable functions
    const range& begin() const { return *this; }
    const range& end() const { return *this; }
    // Iterator functions
    bool operator!=(const range&) const { return iter < last; }
    void operator++() { ++iter; }
    int operator*() const { return iter; }
};

// trial function here
class WF {
  public:

    vector<double> m_psi;   // histogram to be filled
    vector<pos>  m_x;       // walkers
    vector<bool> m_alive;   // walkers

    double m_E0;            // ground state energy
    double m_sumE;          // sum of energy from accepted MC
    double m_sumE2;

    double m_nMCSteps;      // number of total MC steps

    // inputs
    double m_xmin;
    double m_xmax;
    double m_dx;

    double m_dt;              // time step

    int m_N;
    int m_NT;

    int m_reqSteps;

    void clean() {
      m_nMCSteps = 0;
      m_sumE = 0;
      m_sumE2 = 0;
      for (auto &x : m_psi) { x = 0.0; }
    }

    WF(int NT = 10000, int reqSteps = 4000) {
      m_reqSteps = reqSteps;

      // could be inputs
      m_NT = NT;
      m_N = m_NT;

      m_xmin = 0.0;
      m_xmax = 10.0;
      m_dx = 0.1;

      m_dt = 0.05;
      m_E0 = 0;

      m_sumE = 0;
      m_sumE2 = 0;

      m_psi.resize(int((m_xmax - m_xmin)/m_dx));
      m_x.resize(m_N);
      m_alive.resize(m_N);
      uniform_real_distribution<double> rFlat5(-0.5, 0.5);
      for (int i : range(m_N)) {
        m_x[i] = pos();
        for (int d : range(DIM)) {
          m_x[i][d] = rFlat5(rEngine);
        }
        m_alive[i] = true;
      }

      clean();
    }

    double V(pos r) {
      return 0.5*(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));
      //double Z = 0.1;
      //return -Z/sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));
      //if (sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2)) < 2.0) {
      //  return -20.0;
      //}
      //return 0;
    }

    void step(int n) {
      uniform_real_distribution<double> rFlat(0, 1);
      // shift walker
      for (int i : range(DIM)) {
        m_x[n][i] += sqrt(m_dt)*rGaus(rEngine);
      }
      double q = exp(-m_dt*(V(m_x[n]) - m_E0));
      int survivors = int(q);
      if (q - survivors > rFlat(rEngine)) {
        survivors++;
      }
      if (survivors-1 > 0) {
        m_x.resize(m_N+(survivors-1));
        m_alive.resize(m_N+(survivors-1));
        for (int i : range(survivors - 1)) {
          for (int d : range(DIM))
            m_x[m_N][d] = m_x[n][d];

          m_alive[m_N] = true;
          m_N++;
        }
      }
      if (survivors == 0)
        m_alive[n] = false;
      if (m_N > 500000) {
        cout << "Too many nodes: " << m_N << std::endl;
        exit(-1);
      }
    }

    void MC() {
      int N0 = m_N;
      for (int i : range(N0)) step(i);

      int j = 0;
      for (int i : range(m_N)) {
        if (m_alive[i] && i != j) {
          m_x[j] = m_x[i];
          m_alive[j] = m_alive[i];
        }
        if (m_alive[i]) ++j;
      }
      m_N = j;
      m_E0 += log(double(m_NT)/double(m_N))/10.0;
      m_sumE += m_E0;
      m_sumE2 += m_E0*m_E0;
      for (int i : range(m_N)) {
        double r2 = 0;
        for (int d : range(DIM)) {
          r2 += pow(m_x[i][d], 2);
        }
        int k = int(sqrt(r2)/m_xmax*m_psi.size());
        if (k < m_psi.size())
          m_psi[k] += 1;
      }
      m_nMCSteps++;
    }

    void thermalise() {
      for (int i : range(int(0.2*m_reqSteps))) {
        MC();
      }
    }

    void run() {
      thermalise();
      clean();
      for (int i : range(m_reqSteps)) {
        MC();
      }
    }

    double eMean() {
      return m_sumE/double(m_nMCSteps);
    }

    double eError() {
      return sqrt(m_sumE2/double(m_nMCSteps) - pow(m_sumE/double(m_nMCSteps), 2))/sqrt(m_nMCSteps);
    }

    double psiNorm() {
      double n = 0;
      for (int i : range(m_psi.size())) n += pow(i*m_dx, DIM-1)*pow(m_psi[i], 2)*m_dx;
      return n;
    }

    void write(const string &f) {
      ofstream ff(f.c_str());
      ff << "# $E = " << eMean() << " \\pm " << eError() << "$" << endl;
      ff << "# r          psi" << endl;
      double norm = psiNorm();
      for (int i : range(m_psi.size())) {
        //ff << i*m_dx + m_xmin << "    " << pow(i*m_dx, DIM-1)*pow(m_psi[i], 2)/norm << endl;
        ff << i*m_dx + m_xmin << "    " << m_psi[i]/sqrt(norm) << endl;
      }
      ff.close();
    }
};

int main(int argc, char **argv) {

  WF psi(300, 1000);
  psi.run();

  psi.write("psi.dat");

  return 0;
}

