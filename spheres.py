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
    self.bounce = 0
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

'''To take in all values, read lines individually and split, then use a loop to assign
each value to the corresponding component (mass, radius, etc.), where all of these
components are strings or arrays'''

if __name__ == "__main__":

  univ_rad = sys.argv[0]
  univ_coll = sys.argv[1]

  universe = Universe(univ_rad,univ_coll)

  print("Please enter the mass, radius, x/y/z position, x/y/z velocity", "\n", "and name of each sphere", "\n", "When complete, use EOF / Ctrl-D to stop entering")

  name = []
  initials = []
  
  for line in sys.stdin:

    init_val = line.split(" ")

    name.append(init_val[8])

    initials.append(Ball(init_val[0],init_val[1],init_val[8],[init_val[2],init_val[3],init_val[4]],[init_val[5],init_val[6],init_val[7]]))

  # for i in range(len(initials)):
  #   names = initials[i].name.replace("\n", '')
  #   print(names)

''' Below, I write the code established with Taylor setting up the check for
collisions with the universe and all other spheres

The following lines go in the iscollisions statement'''

while ball_array != []:
  for i in range(len(ball_array)):
    if (np.dot(ball_array[i].pos, ball_array[i].pos) >= (radius - ball_array[i].rad)**2):
      Ucollision(ball_array[i])
      ball_array[i].bounce = ball_array[i].bounce + 1

      if (ball_array[i].bounce >= univ_coll):
        print(ball_array[i].names, " has left for the next plane")
        ball_array.remove(i)

      break

    for j in range(i+1,len(ball_array)):
      if (np.dot(np.subtract(ball_array[i].pos,ball_array[j].pos),np.subtract(ball_array[i].pos,ball_array[j].pos)) <= (ball_array[i].rad - ball_array[j].rad)**2):
        update_vel(ball_array[i],ball_array[j])

        if (ball_array[i].bounce >= univ_coll):
          print(ball_array[i].names, " has left for the next plane")
          ball_array.remove(i)

        if (ball_array[j].bounce >= univ_coll):
          print(ball_array[j].names, " has left for the next plane")
          ball_array.remove(j)

        break

''' Below is the mathematical setup for the real collision checks and
time quadratic '''

A = (np.linalg.norm(np.subtract(ball_array[i].pos,ball_array[j].pos)))**2
B = np.dot(np.subtract(ball_array[i].pos,ball_array[j].pos),np.subtract(ball_array[i].vel,ball_array[j].vel))
C = (np.linalg.norm(np.subtract(ball_array[i].pos,ball_array[j].pos)))**2 - (ball_array[i].rad + ball_array[j].rad)**2

# Broken into plus an minus because I couldn't figure out how tf to make +-

t_plus = ((-1) * B + np.sqrt(B ** 2 - A * C)) / A
t_minus = ((-1) * B - np.sqrt(B ** 2 - A * C)) / A