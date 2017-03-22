#include <boost/python.hpp>
#include "QuantumMC.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(quantumMC)
{
  class_<QuantumMC>("QuantumMC", init<int, int>())
    .def("run", &QuantumMC::run)
  ;
}

