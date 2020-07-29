# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 16:44:24 2020
"""
import numpy as np

dt     = .01
N      = 1000  # T/dt
c= 0.7
sig= 0.06

ts    = np.arange(0, N*dt, dt)
ys    = np.zeros(N)
ys[0] = 0
for i in range(1, ts.size):
    t = (i-1) * dt
    dW=np.random.normal(loc=0.0, scale=np.sqrt(dt))
    ys[i] = ys[i-1] - c*ys[i-1]*dt + sig*dW
