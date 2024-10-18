import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from HW6Controller import ctrlPD
from massDynamics1 import massDynamics

mass = massDynamics()
reference = signalGenerator(amplitude=0.5, frequency=0.1)
controller = ctrlPD()
force = signalGenerator(amplitude=10.0, frequency=1)

dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start
while t < P.t_end:
    t_next_plot = t + P.t_plot

    while t < t_next_plot:
        r = reference.square(t)
        x = mass.state
        u = controller.update(r, x)
        y = mass.update(u)
        t = t + P.Ts

    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    t = t + P.t_plot
    plt.pause(0.001)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()