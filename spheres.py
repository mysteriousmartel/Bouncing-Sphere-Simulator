import sys
import numpy

'''For the case of one sphere

print('Please enter initial values on individual lines for each sphere:')

args = sys.stdin.read()

init_val = args.split()

mass = int(init_val[0])
radius = int(init_val[1])
x_0 = int(init_val[2])
y_0 = int(init_val[3])
z_0 = int(init_val[4])
vx_0 = int(init_val[5])
vy_0 = int(init_val[6])
vz_0 = int(init_val[7])

print(mass, " ", radius, " ", x_0, " ", y_0, " ", z_0, " ", vx_0, " ", vy_0, " ", vz_0)'''

'''To take in all values, read lines individually and split, then use a loop to assign
each value to the corresponding component (mass, radius, etc.), where all of these
components are strings or arrays'''

print('Enter for each sphere as: mass, radius, x_0, y_0, z_0, vx_0, vy_0, vz_0')

mass, radius, x_0, y_0, z_0, vx_0, vy_0, vz_0 = # initialize all as strings

for line in sys.stdin:
  init_val = line.split(" ")

  mass = mass.append(init_val[0])
