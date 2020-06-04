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

input_file = "validation_files/input_data.txt"

with open(input_file) as s:
  lines = s.read().splitlines()
  inputs = ("ion_mass_u", "Te", "B0", "R0", "LB", "Ln", "kz", "kx", "N_voltages")
  for line in lines:
    # first I read the input parameters
    # I need to know N_voltages to initialize
    # the numpy arrays.
    temp = line.split()
    try:
      #print(temp[1])
      if temp[0] == 'ion_mass_u':
        ion_mass_u = np.double(temp[1])
      elif temp[0] == "Te":
        Te = np.double(temp[1])
      elif temp[0] == "B0":
        B0 = np.double(temp[1])
      elif temp[0] == "R0":
        R0 = np.double(temp[1])
      elif temp[0] == "plasma_thickness":
        plasma_thickness = np.double(temp[1])
        A_plasma = 2*np.pi*R0*plasma_thickness;
      elif temp[0] == "LB":
        LB = np.double(temp[1])
      elif temp[0] == "Ln":
        Ln = np.double(temp[1])
      elif temp[0] == "kz":
        kz = np.double(temp[1])
      elif temp[0] == "kx":
        kx = np.double(temp[1])
      elif temp[0] == "N_voltages":
        N_voltages = int(temp[1])
      elif temp[0] == "Ni":
        Ni = int(temp[1])
      elif temp[0] == "Nj":
        Nj = int(temp[1])
      elif temp[0] == "tf":
        tf = np.double(temp[1])        
      elif temp[0] == "n10_n0":
        n10_n0 = np.double(temp[1])  
    except(IndexError):
      pass

#################### Validation ######################  
#dir = os.getcwd()  
paths=glob.glob("validation_files/*.txt")
N_v = []
solution_exact = {}
solution_nonlinear = {}
solution_linear = {}
space = {}
err_solver = []

for path in paths:
  if "validation_files/val_space_N" in path:
    N_v.append( int(path[28:-4]) )

for N in N_v:
  templist =[]
  with open("validation_files/val_solution_linear_N{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))

  solution_linear[N] = np.array(templist)
  if N is 100:
    print(solution_linear[N])  
  
  templist =[]
  with open("validation_files/val_solution_nonlinear_N{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))
  solution_nonlinear[N] = np.array(templist)  
  
  templist =[]
  with open("validation_files/val_exact_linear_N{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))
  solution_exact[N] = np.array(templist)   

  templist =[]
  with open("validation_files/val_space_N{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))
  space[N] = np.array(templist)  

# build error dictionary
plt.figure(figsize=(h_size, v_size))

for N in N_v:
  max_index = np.argmax( np.abs (solution_linear[N] - solution_exact[N])  )
  # if N is 100:
    # print(solution_linear[N][0])
  err_solver.append( np.abs( ( solution_linear[N][max_index] - solution_exact[N][max_index]  ) \
                   / solution_exact[N][max_index] ) )
  
N_v = np.array(N_v)
err_solver = np.array(err_solver)

plt.loglog(N_v,err_solver,'o')
plt.title("$n1_0$/$n_0$ = {}".format(n10_n0))
plt.xlabel("Nj")
plt.ylabel("error")
plt.tight_layout() 
plt.savefig("validation_files/Convergence_solver.png") 

for N in N_v:
  plt.figure(figsize=(h_size, v_size))
  plt.plot(space[N],solution_exact[N])
  plt.plot(space[N],solution_linear[N])
  plt.plot(space[N],solution_nonlinear[N])
  plt.title( "N = {}     $n1_0$/$n_0$ = {}".format(N,n10_n0) )
  plt.legend(["exact","linear","nonlinear"])
  plt.xlabel("azimuthal coordinate [m]")
  plt.ylabel("$\phi$ [V]")  
  plt.tight_layout() 
  filename = "validation_files/solution_N{}_{}".format(N,n10_n0)
  filename = filename.replace(".","_") + ".png"
  plt.savefig(filename) 

plt.show()

