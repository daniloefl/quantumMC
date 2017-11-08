
import sys
sys.path.append('lib/')
sys.path.append('../lib/')

import quantumMC
import numpy as np

def V(x):
  return 0.5*x*x

def guidingWF(x):
  #return 1
  return np.exp(-0.5*x*x)

def quantumForce(x): # deriv(guidingWF)/guidingWF
  #return 0
  return -x

def localEnergy(x): # H WF/WF
  #return V(x)
  e = V(x)
  if x != 0: e += -0.5*(-np.exp(-0.5*x*x) + x*x*np.exp(-0.5*x*x))
  return e

c = quantumMC.SchroedingerDiffusionMC(V)
#c.setImportanceSampling(True, guidingWF, localEnergy, quantumForce)
c.setXmin(-5.0)
c.setXmax(5.0)
c.setDeltaX(0.1)
c.setTimeStep(0.001)
c.setNSteps(10000)
c.setN(500)
c.run()
x, psi = c.getPsi()

energy, energy_error = c.getEnergy()
print energy, energy_error

energy, energy_error = c.getLocalEnergy()
print energy, energy_error

from bokeh.plotting import figure, output_file, save
from bokeh.models import Label
output_file("test.html", title="QuantumMC test result")
p = figure(title="QuantumMC test result", x_axis_label='x', y_axis_label='Psi(x)', width = 800, height = 800)

p.circle(x, psi, legend="Psi(x)", line_width=2, line_color="red")
#p.circle(x, [y**2 for y in psi], legend="|Psi(x)|^2", line_width=2, line_color="blue")
energy = Label(x=70, y=700, x_units='screen', y_units='screen',
                 text='Energy = %10.8f +/- %10.8f' % (energy, energy_error), render_mode='css',
                 border_line_color='black', border_line_alpha=0.0,
                 background_fill_color='white', background_fill_alpha=0.0)
p.add_layout(energy)

save(p, "test.html", title = "QuantumMC test result")
