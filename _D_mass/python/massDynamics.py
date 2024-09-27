import numpy as np
class systemDynamics:
    def __init__(self, sample_rate):
        y0 = 0.0
        ydot0 = 0.0
        self.state = np.array([
            [y0],
            [ydot0],
        ])
        self.Ts = sample_rate
        self.limit = 1.0
        self.a0 = 3.0
        self.a1 = 2.0
        self.b0 = 4.0
        alpha = 0.2
        self.a1 = self.a1 * (1.+alpha*(2.*np.random.rand()-1.))
        self.a0 = self.a0 * (1.+alpha*(2.*np.random.rand()-1.))
        self.b0 = self.b0 * (1.+alpha*(2.*np.random.rand()-1.))

    def f(self, state, u):
# for system xdot = f(x,u), return f(x,u)
        y = state.item(0)
        ydot = state.item(1)
# The equations of motion.
        yddot = -self.a1 * ydot - self.a0 * y + self.b0 * u
# build xdot and return
        xdot = np.array([[ydot], [yddot]])
        return xdot

    def h(self):
    # Returns the measured output y = h(x)
        y = self.state.item(0)
 # return output
        return y

    def update(self, u):
# This is the external method that takes the input u(t)
# and returns the output y(t).
        u = self.saturate(u, self.limit) # saturate the input
        self.rk4_step(u) # propagate the state by one time step
        y = self.h() # compute the output at the current state
        return y

    def rk4_step(self, u):
# Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u
