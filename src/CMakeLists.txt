PYTHON_ADD_MODULE(quantumMC QuantumMC.cxx quantumMCPython.cxx)
#add_library(quantumMC SHARED QuantumMC.cxx quantumMCPython.cxx)
target_link_libraries(quantumMC ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

add_executable(dmc_test dmc_test.cxx QuantumMC.cxx)
target_link_libraries(dmc_test ${Boost_LIBRARIES})
