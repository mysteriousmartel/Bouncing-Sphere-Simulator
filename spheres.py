import numpy as np
import sys, os
from decimal import Decimal
import copy
import math
from timeit import default_timer as timer

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
    self.pos += t * self.vel

    for i in range(len(ball_array)):
      ball_array[i].time = 0.00000000

class Universe():

  def __init__(self, radius, max_col, ball_array):
    self.radius=radius
    self.max_col=max_col
    self.ball_array=copy.deepcopy(ball_array)
    self.time = 0.00000000
      
  def __repr__(self):
    str_out = str(self.time) + "\n"
    for i in range(len(self.ball_list)):
      str_out += str(self.ball_list[i].name)+" "+str(self.ball_list[i].pos)[1:-1]+" "+str(self.ball_list[i].vel)[1:-1]
      if i != len(self.ball_list) - 1:
        str_out += "\n"
    return str_out  

  def collideS(self,ball1,ball2):
    delv = np.subtract(ball1.vel,ball2.vel)
    delp = np.subtract(ball1.pos,ball2.pos)
    Radsum = ball1.radius - ball2.radius

    A = np.linalg.norm(delv) ** 2
    B = 2 * np.dot(delp,delv)
    C = np.linalg.norm(delp) ** 2 - Radsum ** 2

    disc = B ** 2 - 4 * A * C

    if (disc < 0 ):
      return None
    else:
      t_plus = ((-1) * B + np.sqrt(disc)) / (2 * A)
      t_minus = ((-1) * B - np.sqrt(disc)) / (2 * A)

      if (t_plus < 0 and t_minus < 0):
        return None
      elif (t_minus < 0):
        return t_plus
      else:
        return t_minus

  def collideU(self,ball,radius):
    A = np.dot(ball.vel,ball.vel)
    B = 2 * np.dot(ball.pos,ball.vel)
    C = np.dot(ball.pos,ball.pos) - radius ** 2

    disc = B ** 2 - 4 * A * C

    if (A == 0):
      return None
    else:
      t_plus = ((-1) * B + np.sqrt(disc)) / (2 * A)
      t_minus = ((-1) * B - np.sqrt(disc)) / (2 * A)

      if (t_plus < 0 and t_minus < 0):
        return None
      elif (t_minus < 0):
        return t_plus
      else:
        return t_minus

  def realColl(self,ball1,ball2,t):
    delv = np.subtract(ball1.vel,ball2.vel)
    r1 = np.sum(ball1.pos,t*np.array(ball1.vel))
    r2 = np.sum(ball2.pos,t*np.array(ball2.vel))
    delr = np.subtract(r1,r2)

    realCheck = np.dot(delr,delv)
    return realCheck
   
  def Ucollision(self, ball1, t):    
    magp=np.linalg.norm(ball1.pos)
    upos=np.array(ball1.pos)*(1/magp)
    
    norm_v = np.dot(ball1.vel,upos)*np.array(upos)
    tan_v = np.subtract(ball1.vel,norm_v)
    ball1.vel = np.subtract(tan_v,norm_v)
    ball1.bounce+=1
    
    print(ball1, " collided with the universe")
    print(ball1,"'s new velocity is ", ball1.vel)
    print("At time ", t)

    for i in range(len(ball_array)):
      ball_array[i].time = t
       
  def Scollision(self, ball1, ball2,t):
    p1, p2= ball1.pos, ball2.pos
    v1, v2= ball1.vel, ball2.vel
    m1, m2 = ball1.mass, ball2.mass
    
    vp=  np.subtract(v1,v2)
    v_p= np.subtract(v2-v1)
    
    rp=  np.subtract(p1-p2)
    r_p= np.subtract(p2-p1)
    
    mp=  np.subtract(m1-m2)
    m_p= np.subtract(m2-m1)
    
    ball1.vel= v1 - (2.*m2/(m1+m2)) * (np.dot(vp,rp)/np.dot(rp,rp))*rp 
    ball2.vel= v2 - (2.*m1/(m1+m2)) * (np.dot(v_p,r_p)/np.dot(r_p,r_p))*r_p 
    
    ball1.bounce+=1
    ball2.bounce+=1
    
    print(ball1, " collided with ", ball2)
    print(ball1,"'s new velocity is ", ball1.vel)
    print(ball2,"'s new velocity is ", ball2.vel)
    print("At time ", t)

    for i in range(len(ball_array)):
      ball_array[i].time = t

  def update_pos(self):
    for ball in self.ball_array:
      ball.updatepos(ball_array[ball].time)
  
  def energy(self, ball_array):
    tot_en=0
    for i in range(len(self.ball_array)):
      tot_en+=.5 * ball_array[i].mass * (np.linalg.norm(ball_array[i].vel))**2
    print("energy: {}".format(tot_en))
  
  def momentum(self, ball_array):
    tot_ro=0
    for i in range(len(self.ball_array)):
      tot_ro+=ball_array[i].mass * np.array(ball_array[i].vel)
    print("momentum: {}".format(tot_ro))

if __name__ == "__main__":

  univ_rad = float(sys.argv[1])
  univ_coll = int(sys.argv[2])

  totalTime = 0.0

  initials = []

  def inputCheck(init_val,mass,radius,x0,y0,z0,vx,vy,vz,name):
    if (len(init_val) != 9):
      return 1
    elif (type(mass) != float and mass <= 0):
      return 1
    elif (type(radius) != float and radius <= 0):
      return 1
    elif (type(x0) != float or type(x0) != int):
      return 1
    elif (type(y0) != float or type(y0) != int):
      return 1
    elif (type(z0) != float or type(z0) != int):
      return 1
    elif (type(vx) != float or type(vx) != int):
      return 1
    elif (type(vy) != float or type(vy) != int):
      return 1
    elif (type(vz) != float or type(vz) != int):
      return 1
    elif (type(name) != str):
      return 1
    else:
      pass

  print("Please enter the mass, radius, x/y/z position, x/y/z velocity", "\n", "and name of each sphere", "\n", "When complete, use EOF / Ctrl-D to stop entering")

  # mass, radius, x0, y0, z0, vx, vy, vz, name = raw_input().split()
  
  for line in sys.stdin:
    init_val = line.split()
    mass = float(init_val[0])
    radius = float(init_val[1])
    x0 = float(init_val[2])
    y0 = float(init_val[3])
    z0 = float(init_val[4])
    vx = float(init_val[5])
    vy = float(init_val[6])
    vz = float(init_val[7])
    name = init_val[8]

    print(type(mass)," ",type(radius)," ",type(x0)," ",type(y0))

    # inputCheck(init_val,mass,radius,x0,y0,z0,vx,vy,vz,name)

    initials.append(Ball(mass,radius,name,[x0,y0,z0],[vx,vy,vz]))

  universe = Universe(univ_rad,univ_coll,initials)
  ball_array = universe.ball_array

# Running the simulation

  while (ball_array != []):
    minty = float("inf")
    colliders = (0,0)

    for i in range(len(ball_array)-1):
      t = universe.collideU(ball_array[i],univ_rad)
      if (t != None and t < minty):
        minty = t
        colliders = (i,-1)
      else:
        for j in range(i+1,len(ball_array)):
          t = universe.collideS(ball_array[i],ball_array[j])
          if (t != None and t < minty):
            collCheck = realColl(ball_array[i],ball_array[j],t)
            if (collCheck < 0):
              minty = t
              colliders = (i,j)

    if (colliders[1] == -1):
      universe.Ucollision(ball_array[colliders[0]], minty)
      universe.energy(ball_array)
      universe.momentum(ball_array)
      ball_array[colliders[0]].bounce += 1

      if (ball_array[colliders[0]].bounce == univ_coll):
        ball_array.remove(colliders[0])
        print(ball_array[colliders[0]].name, " has left")

    else:
      universe.Scollision(ball_array[colliders[0]], ball_array[colliders[1]], minty)
      universe.energy(ball_array)
      universe.momentum(ball_array)
      ball_array[colliders[0]].bounce += 1
      ball_array[colliders[1]].bounce += 1

      if (ball_array[colliders[0]].bounce == univ_coll):
        ball_array.remove(colliders[0])
        print(ball_array[colliders[0]].name, " has left")

      if (ball_array[colliders[1]].bounce == univ_coll):
        ball_array.remove(colliders[1])
        print(ball_array[colliders[1]].name, " has left")

    totalTime += minty

  print("Total time for all spheres to vanish: ", totalTime)