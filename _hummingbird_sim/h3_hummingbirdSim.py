import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
# from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlEquilibrium import ctrlEquilibrium


# instantiate the hummingbird dyanmics
hummingbird = HummingbirdDynamics(alpha=0.0)
ce = ctrlEquilibrium()

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop

    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        ref = np.array([[0.], [0.], [0.]])
        u = ce.update()
        # u = np.array([[pwm_left], [pwm_right]])
        y = hummingbird.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots at rate t_plot
    animation.update(t, hummingbird.state)
    dataPlot.update(t, hummingbird.state, u, ref)

    # the pause causes figure to be displayed during simulation
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
