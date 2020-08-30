import numpy as np
import matplotlib.pyplot as plt
import csv
import glob
import os
import sys

font = {'size'   : 18}
plt.rc('font', **font)

h_size = 9
v_size = 6
linesize = 3.0
linesize_large = 4.0

# subset of modes
unique_m = (3, 4)
do_decr = True
do_incr = False

import numpy as np
import os
import sys



# with open(input_file) as s:
  # lines = s.read().splitlines()
  # inputs = ("ion_mass_u", "Te", "B0", "R0", "LB", "Ln", "kz", "kx", "N_voltages")
  # for line in lines:
    # # first I read the input parameters
    # # I need to know N_voltages to initialize
    # # the numpy arrays.
    # temp = line.split()
    # try:
      # #print(temp[1])
      # if temp[0] == 'ion_mass_u':
        # ion_mass_u = np.double(temp[1])
      # elif temp[0] == "Te":
        # Te = np.double(temp[1])
      # elif temp[0] == "B0":
        # B0 = np.double(temp[1])
      # elif temp[0] == "R0":
        # R0 = np.double(temp[1])
      # elif temp[0] == "plasma_thickness":
        # plasma_thickness = np.double(temp[1])
        # A_plasma = 2*np.pi*R0*plasma_thickness;
      # elif temp[0] == "LB":
        # LB = np.double(temp[1])
      # elif temp[0] == "Ln":
        # Ln = np.double(temp[1])
      # elif temp[0] == "kz":
        # kz = np.double(temp[1])
      # elif temp[0] == "kx":
        # kx = np.double(temp[1])
      # elif temp[0] == "N_voltages":
        # N_voltages = int(temp[1])
      # elif temp[0] == "Ni":
        # Ni = int(temp[1])
      # elif temp[0] == "Nj":
        # Nj = int(temp[1])
      # elif temp[0] == "tf":
        # tf = np.double(temp[1])        
      # elif temp[0] == "n10_n0":
        # n10_n0 = np.double(temp[1])  
    # except(IndexError):
      # pass

n1_file   = "n1.txt"
phi1_file = "phi1.txt"
vx1_file  = "vx1.txt"
vy1_file  = "vy1.txt"

n1 = []
phi1 = []
vx1 = []
vy1 = []

n1_y = []
phi1_y = []

with open(n1_file) as s:
  lines = s.read().splitlines()
  #for line in lines:
  for i in range(len(lines)):
    temp=lines[i].split()
    n1.append(float(temp[0]))
    if i == len(lines)-1:
      for j in range(len(temp)):
        n1_y.append(float(temp[j]))


with open(phi1_file) as s:
  lines = s.read().splitlines()
  #for line in lines:
  for i in range(len(lines)):
    temp=lines[i].split()
    phi1.append(float(temp[0]))
    if i == len(lines)-1:
      for j in range(len(temp)):
        phi1_y.append(float(temp[j]))

with open(vx1_file) as s:
  lines = s.read().splitlines()
  #for line in lines:
  for i in range(len(lines)):
    temp=lines[i].split()
    vx1.append(float(temp[0]))

with open(vy1_file) as s:
  lines = s.read().splitlines()
  #for line in lines:
  for i in range(len(lines)):
    temp=lines[i].split()
    vy1.append(float(temp[0]))

n1 = np.array(n1)
phi1 = np.array(phi1)
vx1 = np.array(vx1)
vy1 = np.array(vy1)
n1_y = np.array(n1_y)
phi1_y = np.array(phi1_y)



t_v = np.linspace(0.0,1e-6,len(n1)) 

plt.figure(figsize=(h_size, v_size))
plt.plot(t_v,n1)

plt.figure(figsize=(h_size, v_size))
plt.plot(t_v,phi1)

plt.figure(figsize=(h_size, v_size))
plt.plot(t_v,vx1)

plt.figure(figsize=(h_size, v_size))
plt.plot(t_v,vy1)

plt.figure(figsize=(h_size, v_size))
plt.plot(n1_y)

plt.figure(figsize=(h_size, v_size))
plt.plot(phi1_y)

plt.show()

