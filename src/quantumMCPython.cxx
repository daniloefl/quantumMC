#include <boost/python.hpp>
#include <string>
#include "SchroedingerDiffusionMC.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(quantumMC)
{
  class_<SchroedingerDiffusionMC>("SchroedingerDiffusionMC", init<object>())
    .def(init<object, double, double, double, int, int>())
    .def("run", &SchroedingerDiffusionMC::run)
    .def("getPsi", &SchroedingerDiffusionMC::getPsi)
    .def("setXmin", &SchroedingerDiffusionMC::setXmin)
    .def("setXmax", &SchroedingerDiffusionMC::setXmax)
    .def("setDeltaX", &SchroedingerDiffusionMC::setDeltaX)
    .def("setNSteps", &SchroedingerDiffusionMC::setNSteps)
    .def("setN", &SchroedingerDiffusionMC::setN)
    .def("getEnergy", &SchroedingerDiffusionMC::getEnergy)
  ;
}

