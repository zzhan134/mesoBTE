import numpy as np
import random
import matplotlib.pyplot as plt
import randomStack as rng
import time

sigma = 0.5
probability_scatter = 0.05

class Particle():
    def __init__(self,xposdist,yposdist,vxdist,vydist):
        self.xpos = xposdist.draw()
        self.ypos = yposdist.draw()
        self.vx = vxdist.draw()
        self.vy = vydist.draw()
    
class Material():
    def __init__(self,xposdistinitial,yposdistinitial,vxdistinitial,vydistinitial):
        self.fig = plt.figure(1)
        plt.ion()
        self.ax = plt.gca() 
        self.ax.set_autoscale_on(False)
        self.ax.axis([0,1,0,1])
        self.fig.show()
        self.e_array = []
        for i in range(0, 5):   #initial_amount of e (user_input)
            a = Particle(xposdistinitial,yposdistinitial,vxdistinitial,vydistinitial)
            self.e_array.append([a, self.ax.plot(a.xpos, a.ypos, 'ro')])

    def electron_amount(self, t):
        y = int(2.5*np.exp(-t))
        return y

    def electron_input(self, t, xposdistinput,yposdistinput,vxdistinput,vydistinput):
        # self.e_array = []
        for i in range(self.electron_amount(t)):
            a = Particle(xposdistinput,yposdistinput,vxdistinput,vydistinput)
            self.e_array.append([a, self.ax.plot(a.xpos, a.ypos, 'ro')])
            
    def scatter(self, dt, reflectdiststack):
        # updating position
        for i in range(len(self.e_array)):
            self.e_array[i][0].xpos = self.e_array[i][0].xpos + self.e_array[i][0].vx * dt
            self.e_array[i][0].ypos = self.e_array[i][0].ypos + self.e_array[i][0].vy * dt
        # bouncing back
            if self.e_array[i][0].ypos >= 1:
                dy = self.e_array[i][0].ypos - 1
                self.e_array[i][0].ypos = self.e_array[i][0].ypos - 2*dy
                self.e_array[i][0].vy = - self.e_array[i][0].vy
            elif self.e_array[i][0].ypos <=0:
                dy = 0 - self.e_array[i][0].ypos
                self.e_array[i][0].ypos = self.e_array[i][0].ypos + 2*dy
                self.e_array[i][0].vy = - self.e_array[i][0].vy
        # update velocity
        length = int(len(self.e_array)* probability_scatter)
        for i in range(length):
            random_index = random.randint(0, len(self.e_array) - 1)
            self.e_array[random_index][0].vx = reflectdiststack.draw()
            self.e_array[random_index][0].vy = reflectdiststack.draw()
        # passing through if x<0 or x>1
        self.e_array = filter(lambda x: (x[0].xpos < 1) and (x[0].xpos > 0), self.e_array)
        print ('\r electrons still simulated: ' + str(len(self.e_array)),)
        for i in range(len(self.e_array)):
            self.e_array[i][1][0].set_xdata(self.e_array[i][0].xpos)
            self.e_array[i][1][0].set_ydata(self.e_array[i][0].ypos)
        
    
    def refreshMat(self,t,dt, xposdistinput,yposdistinput,vxdistinput,vydistinput,reflectdiststack):
        self.electron_input(t,xposdistinput,yposdistinput,vxdistinput,vydistinput)
        self.scatter(dt,reflectdiststack)
  
    def runsim(self,simduration,dt,xposdistinput,yposdistinput,vxdistinput,vydistinput,reflectdiststack):
        t = 0
        while t <= simduration:
            plt.draw()
            t = t + dt
            self.refreshMat(t,dt, xposdistinput,yposdistinput,vxdistinput,vydistinput,reflectdiststack)
            #for i in range(len(mat.e_array)):
            #    a[i].set_ydata(mat.e_array[i].ypos)
            #    a[i].set_xdata(mat.e_array[i].xpos)
            #    b[i].set_xdata([mat.e_array[i].xpos, mat.e_array[i].xpos + mat.e_array[i].vx*dt])
            #    b[i].set_ydata([mat.e_array[i].ypos, mat.e_array[i].ypos + mat.e_array[i].vy*dt])
            time.sleep(dt)

boxdiststack = rng.randomStack(1000,'boxdist.txt')      #square distribution for 0 to .1 for position
boxdist2stack = rng.randomStack(1000,'boxdist2.txt')        #square distribution for 0 to 1 for velocities
reflectdiststack = rng.randomStack(1000,'reflectdist.txt')  #square distribution for -1 to 0 for scattering (unused)

#displaying the distributions from the 'decks'
q = np.zeros(1000)
z = np.zeros(1000)
y = np.zeros(1000)
for i in range(1000):
    q[i] = boxdiststack.draw()
    z[i] = boxdist2stack.draw()
    y[i] = reflectdiststack.draw()
fig2 = plt.figure(2)
h2 = plt.hist(q,np.linspace(-1,1,100),normed=1)
fig2.suptitle('position (x and y pos)')
fig3 = plt.figure(3)
h3 = plt.hist(z,np.linspace(-1,1,100),normed=1)
fig3.suptitle('velocity (vx and vy)')
fig4 = plt.figure(4)
h4 = plt.hist(y,np.linspace(-1,1,100),normed=1)
fig4.suptitle('scattering (vx and vy)')
fig2.show()
fig3.show()
fig4.show()
u = input("enter any number to continue to the simulation")
plt.close('all')

#simulation time! notice all electrons start in bottom left, and go towards upper right of window as predicted by distributions

mat = Material(boxdiststack,boxdiststack,boxdist2stack,boxdist2stack) #when initializing the material, 
#electrons initialized according to distributions given in argument
#argument format (xpos,ypos,vx,vy)
mat.runsim(20,0.01,boxdiststack,boxdiststack,boxdist2stack,boxdist2stack,reflectdiststack) #when running sim, 
#electrons are added according to distributions given in argument
#argument format (simduration,dt,xpos,ypos,vx,vy)
'''
points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
from scipy.spatial import Delaunay
tri = Delaunay(points)
plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.show()
'''