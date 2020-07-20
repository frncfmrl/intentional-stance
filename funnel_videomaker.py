# -*- coding: utf-8 -*-

"""
Created on Mon Jul 20 12:30:34 2020
 
"""
 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from utilities import sign
         
class funnel_simulation:

    def __init__(self,num_of_dots):
         
         self.nDots=num_of_dots
         
         self.jitterRatio_x=0 #but differerent dots might jitter at a different noise level
         self.jitterRatio_y=0

         self.frame = {
            "yMin": -1.5,
            "yMax": 1.5,
            "xMin": -1,
            "xMax": +1,
            }
         
         self.bounds=(self.frame['xMin'],self.frame['xMin'],self.frame['yMin'],self.frame['yMax'])
         self.Lx=self.frame ['xMax']-self.frame ['xMin']
         self.Ly=self.frame ['yMax']-self.frame ['yMin']
         self.sx=self.jitterRatio_x*self.Lx
         self.sy=self.jitterRatio_y*self.Ly
         
         self.init_state=np.empty((nDots,2))  #initialize
         self.init_state[:,0]=self.frame['xMin']+self.Lx*np.random.rand(nDots) 
         self.init_state[:,1]=self.frame['yMin']+self.Ly*np.random.rand(nDots)
         for n in range(nDots):
              self.init_state[n,:]=self.keep_within(self.init_state[n,:])

         self.state = self.init_state.copy()
         self.time_elapsed = 0
         #tMax=3
         #self.r=self.trajectory(nDots=nDots,dt=dt,nSteps=int(np.floor(tMax/dt)),sx=sx,sy=sy)

    def x_of_y(self,y):
          x= np.exp(y-1)
          return(x)
              
    def keep_within(self,r):
          assert(len(r)==2)
          if abs(r[0])>self.x_of_y(r[1]):
               r[0] =.99*sign(r[0])*self.x_of_y(r[1])
          return(r)
                        
    def draw_manifold(self):
          ys=np.arange(self.frame["yMin"], self.frame["yMax"],.05*(self.frame["yMax"]-self.frame["yMin"]))
          plt.plot(self.x_of_y(ys),ys, color='black', linestyle='solid')
          plt.plot(-self.x_of_y(ys),ys, color='black', linestyle='solid')
          return()
                             
    def repulsion(self,d):
          assert(d>=0) # argument is a distance so it must be positive
          v=1/(d+1)
          return(v)
                                  
    def evolve(self,r,dt,x_of_y,sig_x,sig_y):
          vx=-sign(r[0])*self.repulsion(self.x_of_y(r[1])-abs(r[0]))
          vy=-5
          r[0]=r[0]+vx*dt+sig_x*np.random.randn(1)
          r[1]=r[1]+vy*dt+sig_y*np.random.rand(1)
          return(r)
                                       
    def trajectory(self,nDots=1,dt=0.1,nSteps=100,sx=0,sy=0):
                                            
          # initializing the positions array
          r=np.empty((2,nDots,nSteps)) 
          r[:,:,0]=self.init_state
          for n in range(nDots):
               r[n,:,0]=self.keep_within(r[n,:0])
               for tn in range(nSteps-1):
                    # print(tn)
                    r[n,:,tn+1]=self.evolve(r[n,:,tn],dt,self.x_of_y,sx,sy)
                    r[n,:,tn+1]=self.keep_within(r[n,:,tn+1])
                    #assert(abs(r[n,0,tn+1])<self.x_of_y([n,1,tn+1]))
                    
          return(r)
           
    def showManifold(self): 
              
         # ys=np.arange(frame["yMin"], frame["yMax"],.05*Ly)
         # rightBorder = ax.plot(self.x_of_y(ys),ys, color='black',linestyle='solid')
         # leftBorder = ax.plot(-self.x_of_y(ys),ys, color='black',linestyle='solid')
         #plt.plot(r[0,0,:],r[1,0,:])
         return()
      
    def step(self, dt):
       
         """step once by dt seconds"""
         self.time_elapsed += dt

         for n in range(self.nDots):
               self.state[n,:]=self.evolve(self.state[n,:],dt,self.x_of_y,self.sx,self.sy)
               self.state[n,:]=self.keep_within(self.state[n,:])
          #assert(abs(r[0,n,tn+1])<self.x_of_y([1,n,tn+1]))self.state[:, :2] += dt * self.state[:, 2:]


#------------------------------------------------------------

##################### SET UP SIMULATION  ######################

nDots=3
dt = 1/30 # 30fps
simulation = funnel_simulation(nDots)


######################## CREATE MOVIE OF TRAJECTORIES #######################

fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(simulation.frame['xMin'], simulation.frame['xMax']), 
                     ylim=(simulation.frame['yMin'], simulation.frame['yMax']))

particles, = ax.plot([], [], 'bo', ms=6) # particles holds the locations of the particles
rect = plt.Rectangle(simulation.bounds[::2],
                     simulation.bounds[1] - simulation.bounds[0],
                     simulation.bounds[3] - simulation.bounds[2],
                     ec='none', lw=2, fc='none') # rect is the box edge
ax.add_patch(rect)

def init():
    """initialize animation"""
    global simulation, rect
    particles.set_data([], [])
    rect.set_edgecolor('none')
    return(particles, rect)

def animate(i):
    """perform animation step"""
    global simulation, rect, dt, ax, fig
    simulation.step(dt)

## update pieces of the animation
    rect.set_edgecolor('k')
    particles.set_data(simulation.state[0,:], simulation.state[1,:])
    ms=8 # ms = int(fig.dpi * 2 * simulation.size * fig.get_figwidth() / np.diff(ax.get_xbound())[0])
    particles.set_markersize(ms)
    return(particles, rect)

nFrames=1000
ani = animation.FuncAnimation(fig, animate, frames=nFrames, interval=10, blit=True, init_func=init)
     
# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

ani.save('../funnel_video.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()

################## CREATE STATIC IMAGE OF TRAJECTORIES ############### 

simulation2 = funnel_simulation(nDots)
trajectories=np.empty((nDots,2,nFrames)) 
     
for frameInd in range(nFrames):
     trajectories[:,:,frameInd]=simulation2.state
     #animate(frameInd)
     simulation2.step(dt)
     
# plt.plot(trajectories[0,0,:])
     
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(simulation.frame['xMin'], simulation.frame['xMax']), 
                     ylim=(simulation.frame['yMin'], simulation.frame['yMax']))

rect = plt.Rectangle(simulation.bounds[::2],
                     simulation.bounds[1] - simulation.bounds[0],
                     simulation.bounds[3] - simulation.bounds[2],
                     ec='none', lw=2, fc='none') # rect is the box edge
ax.add_patch(rect)

for n in range(nDots):
     plt.scatter(trajectories[n,0,:], trajectories[n,1,:])

# also plot funnel:
ys=np.arange(simulation.frame["yMin"], simulation.frame["yMax"],.05*(simulation.frame["yMax"]-simulation.frame["yMin"]))
plt.plot(simulation.x_of_y(ys),ys, color='black', linestyle='solid')
plt.plot(-simulation.x_of_y(ys),ys, color='black', linestyle='solid')
 
