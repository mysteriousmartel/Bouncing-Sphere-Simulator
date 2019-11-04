import numpy as np
import sys, os
from decimal import Decimal
import copy
import math
from timeit import default_timer as timer


#creates class for making the balls


class Ball():
  def __init__(self, mass, radius, name, pos=[], vel=[]):
    self.mass = mass
    self.radius = radius
    self.name = name
    self.pos = np.array(pos, np.float64)
    self.vel = np.array(vel, np.float64)
    self.time = 0.00000000

  def __repr__(self):
    return str(self.name)+" "+str(self.pos)[1:-1]+" "+str(self.vel)[1:-1]

  def updatePos(self, t): 
    # dt = (t - self.time)
    dt = float(Decimal(t) - Decimal(self.time))
    self.pos += dt * self.vel
    self.time = t

    
    # Universe Class
    
      #constructor
      #collision checker
      #collision time
      #collisions 
      #update ball position
      #update ball velocity
      #update time
      
    
    
