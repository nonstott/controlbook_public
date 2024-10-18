import numpy as np
import massParam as P

class ctrlPD: 
    def __init__(self):
        # PD gains
        self.kp = 12
        self.kd = 4.5
        #print(’kp: ’, self.kp)
        #print(’kd: ’, self.kd)

    def update(self, z_r, x):
        z = x[0][0]
        zdot = x[1][0]
      
        F_tilde = self.kp * (z_r - z) - self.kd * zdot
         # compute total torque
        F = F_tilde
         #tau = tau_e + tau_tilde
         # always saturate to protect hardware
        F = self.saturate(F, P.F_max)
        return F

    def saturate(u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u