# -*- coding: utf-8 -*-

"""
Created on Mon Jul 20 12:30:34 2020 
Based on :  https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.animation.FuncAnimation.html

 Saves the animation as an mp4.  This requires ffmpeg or mencoder to be
 installed.  The extra_args ensure that the x264 codec is used, so that
 the video can be embedded in html5.  You may need to adjust this for
 your system: for more information, see

 http://matplotlib.sourceforge.net/api/animation_api.html
     
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from utilities import sign
from datetime import datetime

#

class funnel_simulation:

    def __init__(self,num_of_dots,theta):
         
         assert(theta<360)
         self.nDots=num_of_dots
         self.theta=theta*2*np.pi/360
         
         self.jitterRatio_x=0 #but differerent dots might jitter at a different noise level
         self.jitterRatio_y=0

         self.frame = {
            "yMin": -1, #-.5
            "yMax": 1, # 1.5
            "xMin": -1,
            "xMax": +1,
            }
         
         self.bounds=(self.frame['xMin'],self.frame['xMin'],self.frame['yMin'],self.frame['yMax'])
         self.Lx=self.frame ['xMax']-self.frame ['xMin']
         self.Ly=self.frame ['yMax']-self.frame ['yMin']
         self.sx=self.jitterRatio_x*self.Lx
         self.sy=self.jitterRatio_y*self.Ly
         
         self.init_state=np.empty((self.nDots,2))  #initialize
         self.init_state[:,0]=self.frame['xMin']+self.Lx*np.random.rand(nDots) 
         self.init_state[:,1]=self.frame['yMin']+self.Ly*np.random.rand(nDots)
         
         for n in range(nDots):
              self.init_state[n,:]=self.keep_within(self.init_state[n,:])
              self.init_state[n,:]=self.rotate(self.init_state[n,:])
                            
         self.state = self.init_state.copy()
         self.time_elapsed = 0

    def rotate(self,r):
         R=np.zeros((2,2))
         R[0,0]=np.cos(self.theta) 
         R[0,1]=-np.sin(self.theta)
         R[1,0]=np.sin(self.theta)
         R[1,1]=np.cos(self.theta)
         return(R@r)

    def antirotate(self,r):
         R=np.zeros((2,2))
         R[0,0]=np.cos(self.theta) 
         R[0,1]=np.sin(self.theta)
         R[1,0]=-np.sin(self.theta)
         R[1,1]=np.cos(self.theta)
         return(R@r)
         
    def x_of_y(self,y):
          x= np.exp(y-0.5)
          return(x)
              
    def keep_within(self,r):
          assert(len(r)==2)
          if abs(r[0])>self.x_of_y(r[1]):
               r[0] =.95*sign(r[0])*self.x_of_y(r[1])
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

    def evolve(self,r,dt,x_of_y):
          vx=-sign(r[0])*self.repulsion(self.x_of_y(r[1])-abs(r[0]))
          vy=-5+2*np.random.randn(1)
          r[0]=r[0]+vx*dt
          r[1]=r[1]+vy*dt
          return(r)

    def trajectory(self,nDots=1,dt=0.1,nSteps=100,sig_x=0,sig_y=0):
          # initializing the positions array
          r=np.empty((2,nDots,nSteps)) 
          r[:,:,0]=self.init_state                    
          for n in range(nDots):
               r[n,:,0]=self.keep_within(r[n,:0])
               for tn in range(nSteps-1):
                    r[n,:,tn+1]=self.antirotate(r[n,:,tn])
                    r[n,:,tn+1]=self.evolve(r[n,:,tn+1],dt,self.x_of_y)
                    r[n,:,tn+1]=self.keep_within(r[n,:,tn+1])
                    r[n,:,tn+1]=self.rotate(r[n,:,tn+1])
                    #assert(abs(r[n,0,tn+1])<self.x_of_y([n,1,tn+1]))
          return(r)
          
    def manifold(self): 
         #
         # ys=np.arange(self.frame["yMin"], self.frame["yMax"],.05*self.Ly)
         #
         ys=np.arange(-np.sqrt(2),np.sqrt(2),.05)
         lWall=np.array([simulation.x_of_y(ys), ys])
         rWall=np.array([-simulation.x_of_y(ys), ys])
         lWall=self.rotate(lWall)
         rWall=self.rotate(rWall)
         #
         return(rWall[0],rWall[1],lWall[0],lWall[1])
      
    def step(self,dt,deltas):
         """step once by dt seconds"""
         self.time_elapsed += dt
         for n in range(self.nDots):
               self.state[n,:]=self.antirotate(self.state[n,:])
               self.state[n,:]=self.evolve(self.state[n,:],dt,self.x_of_y)+deltas[n,:]
               self.state[n,:]=self.keep_within(self.state[n,:])
               self.state[n,:]=self.rotate(self.state[n,:])
               #assert(abs(r[n,0,tn+1])<self.x_of_y([n,1,tn+1]))

############################### MODULE EXECUTION ########################################

if __name__=="__main__":  
     
     ######################### SIMULATION PARAMETERS #########################
     
     ## set up simulation 
     nDots=10
     nFrames=450
     dt = 1/750 # 30fps
     angle=np.random.rand(1)[0]*360
     
     ## set up movie
     simulation = funnel_simulation(nDots,angle)
     fig = plt.figure()
     fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
     
     ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                          xlim=(simulation.frame['xMin'], simulation.frame['xMax']), 
                          ylim=(simulation.frame['yMin'], simulation.frame['yMax']))
     
     ax.axis('off')
     ax.set_xlim(simulation.frame['xMin'], simulation.frame['xMax'])
     ax.set_ylim(simulation.frame['yMin'], simulation.frame['yMax'])
     
     rightWall, =ax.plot([],[], color='black', linestyle='solid')
     leftWall, =ax.plot([],[], color='black', linestyle='solid')
     particles, = ax.plot([], [], 'bo', ms=6) # particles holds the locations of the particles
     
     def init():
         #
         """initialize animation"""
         #
         global simulation, particles,rightWall,leftWall
         rightWall.set_data([],[])    
         rightWall.set_data([],[])    
         particles.set_data([], [])
         #
         return particles,rightWall,leftWall,
     
     def animate(i):
     
         """perform animation step"""
         
         global simulation,particles,rightWall,leftWall,dt,ax,fig,rightWall_xs,rightWall_ys,nFrames,all_deltas 
         simulation.step(dt,deltas=all_deltas[i,:,:]) #all_deltas must have size (nFrames,nDots,2)
     
         #update dots
         particles.set_data(simulation.state[:,0],simulation.state[:,1])
         ms=8 # ms = int(fig.dpi * 2 * simulation.size * fig.get_figwidth() / np.diff(ax.get_xbound())[0])
         particles.set_markersize(ms)
         
         # update walls

         unrotated_yMax=-100
         for n in range(simulation.nDots):
              ss=simulation.antirotate(simulation.state[n,:])
              unrotated_yMax=np.max((unrotated_yMax,ss[1]))
              
         #if unrotated_yMax<(simulation.frame["yMin"]+simulation.frame["yMax"])/2:
         if i>nFrames/2:
              rWall_x,rWall_y,lWall_x,lWall_y=simulation.manifold()
              rightWall.set_data(rWall_x,rWall_y)
              leftWall.set_data(lWall_x,lWall_y)
              
         return particles,rightWall,leftWall,
     
     def ornstUhlen(nTimes,dt=0.05,c=0.7,sig=0.06):
          #
          #outputs is an arry of length (nTimes)
          #
          ts    = np.arange(0, nTimes*dt, dt)
          ys    = np.zeros(nTimes)
          ys[0] = 0
          for i in range(1, ts.size):
               dW=np.random.normal(loc=0.0, scale=np.sqrt(dt))
               ys[i] = ys[i-1] - c*ys[i-1]*dt + sig*dW
          #
          return(ys)
     
     ######################### MAKE AND SAVE VIDEOS #########################
     
     #create fluctuation matrix of size (nFrames,nDots,2) 
     # To add: Let's compute the numeber of Frames as the number of iteration after which all the dots are below the figure's edge
     
     all_deltas=np.empty((nFrames,nDots,2))
     for n in range(nDots):
          all_deltas[:,n,0]=ornstUhlen(nFrames,dt=dt,c=0.7,sig=.5)
          all_deltas[:,n,1]=ornstUhlen(nFrames,dt=dt,c=0.7,sig=0)
     
     ani = animation.FuncAnimation(fig,animate,np.arange(1, nFrames),init_func=init,interval=25,blit=True)
     
     # ani = animation.FuncAnimation(fig,animate,frames=nFrames,interval=10,blit=True,init_func=init)
     # If blit == True, the "animate" function must return an iterable of all artists that were modified or created. 
     
     timeNow=datetime.now()
     timestr=timeNow.strftime("%y%m%d_%H%M%S")
     videoName= "../outputs/funnelVideo_"+timestr+".mp4"
     
     ani.save(videoName, fps=30, extra_args=['-vcodec', 'libx264'])
     
# might add the functionality to make the external border of the movie circular. 