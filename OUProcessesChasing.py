# Simple OU process, notation from http://web.math.ku.dk/~susanne/StatDiff/Overheads1b

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

T = 1
dt = 0.005
iterations = int(T//dt)

realisations = 20

beta = 5.
alpha = 0.
sigma = 50.

dW = sigma * np.random.randn(realisations, iterations)
x = 50. * np.random.randn(realisations, iterations)

wall_x = np.linspace(0, T, iterations)
wall1_y = np.exp(np.log(np.max(x[:,0])) - beta*wall_x + alpha)
wall2_y = -np.exp(np.log(np.max(x[:,0])) - beta*wall_x + alpha)

initial_spread = 7.
initial_time = np.random.randn(realisations, 1)*initial_spread

fig = plt.figure(figsize=(10,15))
ax1 = fig.add_subplot(111, autoscale_on=True, ylim=(-T, 5*dt), xlim=(-np.max(x[:,0])*2, np.max(x[:,0])*2))
ax1.grid()
time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

particles, = ax1.plot([], [], 'o', lw=2, label='Particle')

wall1, = ax1.plot([], [], '-k', lw=2)
wall2, = ax1.plot([], [], '-k', lw=2)

def init():
    """initialize animation"""
    time_text.set_text('')
    particles.set_data([], [])

    wall1.set_data([], [])
    wall2.set_data([], [])

    return particles, time_text, wall1, wall2

def animate(i, wall1):
    dx = - beta * (x[:,i] - alpha) + dW[:,i]/np.sqrt(dt)
    x[:,i+1] = x[:,i] + dt * dx

    time_text.set_text('iteration = %.2f' % i)
    particles.set_data([x[:,i]], [-(initial_time+i)*dt])

    if i > iterations/4:
        wall1.set_data([wall1_y], [-wall_x])
        wall2.set_data([wall2_y], [-wall_x])

    return particles, time_text, wall1, wall2

ani = animation.FuncAnimation(fig, animate, fargs=[wall1], frames=iterations,
                              interval=1, blit=True, init_func=init)

# Set up formatting for the movie files
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

# ani.save('test.mp4', writer=writer)


plt.show()