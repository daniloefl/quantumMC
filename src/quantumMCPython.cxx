#include <boost/python.hpp>
#include "QuantumMC.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(quantumMC)
{
  class_<QuantumMC>("QuantumMC", init<>())
    .def(init<double, double, double, int, int>())
    .def("run", &QuantumMC::run)
    .def("getPsi", &QuantumMC::getPsi)
    .def("setXmin", &QuantumMC::setXmin)
    .def("setXmax", &QuantumMC::setXmax)
    .def("setDeltaX", &QuantumMC::setDeltaX)
    .def("setNSteps", &QuantumMC::setNSteps)
    .def("setN", &QuantumMC::setN)
  ;
}

