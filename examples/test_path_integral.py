import sys
sys.path.append('lib/')
sys.path.append('../lib/')

def V(x):
  return 0.5*x*x

def dVdx(x):
  return x

import quantumMC
c = quantumMC.PathIntegralMC(V, dVdx)
c.setNSteps(10000)
c.run()
x, psi = c.getPsi()
energy, energy_error = c.getEnergy()

from bokeh.plotting import figure, output_file, save
from bokeh.models import Label
output_file("test_path_integral.html", title="Path integral simulation test result")
p = figure(title="Path integral test result with harmonic oscillator", x_axis_label='x', y_axis_label='Psi(x)', width = 800, height = 800)

p.circle(x, psi, legend="Psi(x)", line_width=2, line_color="red")
p.circle(x, [y**2 for y in psi], legend="|Psi(x)|^2", line_width=2, line_color="blue")
energy = Label(x=70, y=750, x_units='screen', y_units='screen',
                 text='Energy = %.5f +/- %.5f' % (energy, energy_error), render_mode='css',
                 border_line_color='black', border_line_alpha=0.0,
                 background_fill_color='white', background_fill_alpha=0.0)
p.add_layout(energy)

save(p, "test_path_integral.html", title = "Path integral test result with harmonic oscillator")
