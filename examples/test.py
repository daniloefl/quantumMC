
import quantumMC
c = quantumMC.QuantumMC(300, 4000)
c.run()
x, psi = c.getPsi()

from bokeh.plotting import figure, output_file, show
output_file("test.html", title="QuantumMC test result")
p = figure(title="QuantumMC test result", x_axis_label='x', y_axis_label='Psi(x)')

p.line(x, psi, legend="Psi(x)", line_width=2, line_color="red")
p.line(x, psi**2, legend="|Psi(x)|^2", line_width=2, line_color="blue")

save(p, "test.html", title = "QuantumMC test result")
