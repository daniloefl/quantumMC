PYTHON_ADD_MODULE(quantumMC QuantumMC.cxx quantumMCPython.cxx)
target_link_libraries(quantumMC ${Boost_LIBRARIES} boost_python)

