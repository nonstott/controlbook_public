import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter

reference = signalGenerator(amplitude=0.5, frequency=0.1)
zSig = signalGenerator(amplitude=0.05, frequency=0.5, y_offset=0.2)
thetaSig = signalGenerator(amplitude=np.pi/8, frequency=0.1, y_offset = 0.0)
fSig = signalGenerator(amplitude=5, frequency=.5)

dataPlot = dataPlotter()
animation = blockbeamAnimation()

t= P.t_start
while t < P.t_end:
    r = reference.square(t)
    z = zSig.sin(t)
    theta = thetaSig.sin(t)
    f = fSig.sawtooth(t)

    state = np.array([[0.25+z], [theta], [0.0], [0.0]])
    animation.update(state)
    dataPlot.update(t, r, state, f)

    t = t + P.t_plot
    plt.pause(0.1)

print('Press key to close')
plt.waitforbuttonpress()
plt.close()
