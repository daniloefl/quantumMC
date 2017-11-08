#include <boost/python.hpp>
#include <string>
#include "SchroedingerDiffusionMC.h"
#include "PathIntegralMC.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(quantumMC)
{
  class_<SchroedingerDiffusionMC>("SchroedingerDiffusionMC", init<object>())
    .def(init<object, double, double, double, int, int>())
    .def("setImportanceSampling", &SchroedingerDiffusionMC::setImportanceSampling)
    .def("run", &SchroedingerDiffusionMC::run)
    .def("getPsi", &SchroedingerDiffusionMC::getPsi)
    .def("setTimeStep", &SchroedingerDiffusionMC::setTimeStep)
    .def("setXmin", &SchroedingerDiffusionMC::setXmin)
    .def("setXmax", &SchroedingerDiffusionMC::setXmax)
    .def("setDeltaX", &SchroedingerDiffusionMC::setDeltaX)
    .def("setNSteps", &SchroedingerDiffusionMC::setNSteps)
    .def("setN", &SchroedingerDiffusionMC::setN)
    .def("logGrid", &SchroedingerDiffusionMC::logGrid)
    .def("getEnergy", &SchroedingerDiffusionMC::getEnergy)
    .def("getLocalEnergy", &SchroedingerDiffusionMC::getLocalEnergy)
  ;

  class_<PathIntegralMC>("PathIntegralMC", init<object, object>())
    .def(init<object, object, double, double, double, double, int, double, int>())
    .def("run", &PathIntegralMC::run)
    .def("getPsi", &PathIntegralMC::getPsi)
    .def("setXmin", &PathIntegralMC::setXmin)
    .def("setXmax", &PathIntegralMC::setXmax)
    .def("setDeltaX", &PathIntegralMC::setDeltaX)
    .def("setDelta", &PathIntegralMC::setDelta)
    .def("setDeltaT", &PathIntegralMC::setDeltaT)
    .def("setNSteps", &PathIntegralMC::setNSteps)
    .def("getEnergy", &PathIntegralMC::getEnergy)
  ;
}

