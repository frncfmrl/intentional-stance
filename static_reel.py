# -*- coding: utf-8 -*-

"""
Created on Mon Jul 27 20:30:36 2020
"""

################## CREATE STATIC IMAGE OF TRAJECTORIES ############### 



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from datetime import datetime
from funnel_videomaker import funnel_simulation
     
simulation= funnel_simulation(nDots,0)
trajectories=np.empty((nDots,2,nFrames)) 
     
for frameInd in range(nFrames):
     trajectories[:,:,frameInd]=simulation.state
     #animate(frameInd)
     simulation.step(dt)
     
# plt.plot(trajectories[0,0,:])
     
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(simulation.frame['xMin'], simulation.frame['xMax']), 
                     ylim=(simulation.frame['yMin'], simulation.frame['yMax']))
ax.axis('off')

for n in range(nDots):
     plt.scatter(trajectories[n,0,:], trajectories[n,1,:])

ys=np.arange(simulation.frame["yMin"], simulation.frame["yMax"],\
             .05*(simulation.frame["yMax"]-simulation.frame["yMin"]))

ax.plot(simulation.x_of_y(ys),ys, color='black', linestyle='solid')
ax.plot(-simulation.x_of_y(ys),ys, color='black', linestyle='solid')

ax.set_xlim(simulation.frame['xMin'], simulation.frame['xMax'])
ax.set_ylim(simulation.frame['yMin'], simulation.frame['yMax'])


timeNow=datetime.now()
timestr=timeNow.strftime("%y%m%d_%H%M%S")

plt.show()

figureName= "../outputs/trajectories_"+timestr+".png"
plt.savefig(figureName, bbox_inches='tight')
