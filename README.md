# quantumMC

Small project to illustrate how a diffuse MCMC
can be made to solve the Schroedinger equation.

It creates an interface with Python and allows one to use Python
to steer the MCMC C++ code.
Take a look at examples/test.py, for instance.

# Compilation

This needs boost libraries range and python. Install them with:

```
sudo apt install libboost-dev* libboost-python*
```

To compile:

```
cmake .
make
```


