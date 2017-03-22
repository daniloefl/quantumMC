#include "QuantumMC.h"

int main(int argc, char **argv) {
  QuantumMC qmc(300, 1000);
  qmc.run();
  qmc.write("result.txt");
  return 0;
}
