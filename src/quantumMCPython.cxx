#include <boost/python.hpp>
#include <string>
#include "QuantumMC.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(quantumMC)
{
  class_<QuantumMC>("QuantumMC", init<std::string>())
    .def(init<std::string, double, double, double, int, int>())
    .def("run", &QuantumMC::run)
    .def("getPsi", &QuantumMC::getPsi)
    .def("setXmin", &QuantumMC::setXmin)
    .def("setXmax", &QuantumMC::setXmax)
    .def("setDeltaX", &QuantumMC::setDeltaX)
    .def("setNSteps", &QuantumMC::setNSteps)
    .def("setN", &QuantumMC::setN)
    .def("getEnergy", &QuantumMC::getEnergy)
  ;
}

