import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter

reference = signalGenerator(amplitude=0.5, frequency=0.1)
zSig = signalGenerator(amplitude=1, frequency=0.5, y_offset=0.2)
fSig = signalGenerator(amplitude=2, frequency=0.5)

dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start
while t < P.t_end:
    r = reference.square(t)
    z = zSig.sin(t)
    f = fSig.sawtooth(t)
    state = np.array([[z], [0.0]])
    animation.update(state)
    dataPlot.update(t, r, state, f)

    t = t + P.t_plot
    plt.pause(0.1)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()