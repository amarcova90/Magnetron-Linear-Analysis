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


##### modified wavenumber analysis ######
def kpdy(kdy):
  return (-2.0*np.sin(2*kdy) + 16.0*np.sin(kdy))/12.0

Np = 1000
kdy_v = np.linspace(0.0, np.pi, Np)

plt.figure()
plt.plot(kdy_v,kdy_v)
plt.plot(kdy_v, kpdy(kdy_v))
plt.legend(["exact","$4^{th}$ central"])
plt.xlabel("$k \Delta y$")
plt.ylabel("$k' \Delta y$")
plt.tight_layout() 

max = np.max(kpdy(kdy_v))
print("maximum k'dy = {}".format(max))

input_file = "validation_files/input_data.txt"

with open(input_file) as s:
  lines = s.read().splitlines()
  inputs = ("ion_mass_u", "Te", "B0", "R0", "LB", "Ln", "kz", "kx", "Nj_voltages")
  for line in lines:
    # first I read the input parameters
    # I need to know Nj_voltages to initialize
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
      elif temp[0] == "Nj_voltages":
        Nj_voltages = int(temp[1])
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
Nj_v = []
solution_exact = {}
solution_nonlinear = {}
solution_linear = {}
space = {}
err_solver = []


for path in paths:
  if "validation_files/val_space_Nj" in path:
    Nj_v.append( int(path[29:-4]) )

for N in Nj_v:
  templist =[]
  with open("validation_files/val_solution_linear_Nj{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))

  solution_linear[N] = np.array(templist)
  # if N is 100:
    # print(solution_linear[N])  
  
  templist =[]
  with open("validation_files/val_solution_nonlinear_Nj{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))
  solution_nonlinear[N] = np.array(templist)  
  
  templist =[]
  with open("validation_files/val_exact_linear_Nj{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))
  solution_exact[N] = np.array(templist)   

  templist =[]
  with open("validation_files/val_space_Nj{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))
  space[N] = np.array(templist)  

  plt.figure(figsize=(h_size, v_size))
  plt.plot(space[N],solution_exact[N])
  plt.plot(space[N],solution_linear[N])
  plt.plot(space[N],solution_nonlinear[N])
  plt.title( "$N_j$ = {}     $n1_0$/$n_0$ = {}".format(N,n10_n0) )
  plt.legend(["exact","linear","nonlinear"])
  plt.xlabel("azimuthal coordinate [m]")
  plt.ylabel("$\phi$ [V]")  
  plt.tight_layout() 
  filename = "validation_files/solution_Nj{}_{}".format(N,n10_n0)
  filename = filename.replace(".","_") + ".png"
  plt.savefig(filename) 

for N in Nj_v:
  max_index = np.argmax( np.abs (solution_linear[N] - solution_exact[N])  )
  # if N is 100:
    # print(solution_linear[N][0])
  err_solver.append( np.abs( ( solution_linear[N][max_index] - solution_exact[N][max_index]  ) \
                   / solution_exact[N][max_index] ) )

Nj_v = np.array(Nj_v)
err_solver = np.array(err_solver)

plt.figure(figsize=(h_size, v_size))
plt.loglog(Nj_v,err_solver,'o')
plt.title("$n1_0$/$n_0$ = {}".format(n10_n0))
plt.xlabel("$N_j$")
plt.ylabel("error")
plt.tight_layout() 
plt.savefig("validation_files/Convergence_solver.png") 

# RK4 validation

Ni_v = []
solution_exact_RK4 = {}
solution_nonlinear_RK4 = {}
solution_linear_RK4 = {}
time_RK4 = {}
err_RK4 = []


# timeRK4 = []
# solution_exact_RK4 = []
# solution_linear_RK4 = []
# solution_nonlinear_RK4 = []

for path in paths:
  if "validation_files/valRK4_time_Ni" in path:
    Ni_v.append( int(path[31:-4]) )


for N in Ni_v:
  templist = []
  with open("validation_files/valRK4_exact_linear_Ni{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))

  solution_exact_RK4[N] = np.array(templist)

  templist = []
  with open("validation_files/valRK4_solution_linear_Ni{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))

  solution_linear_RK4[N] = np.array(templist)

  templist = []
  with open("validation_files/valRK4_solution_nonlinear_Ni{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))

  solution_nonlinear_RK4[N] = np.array(templist)

  templist = []
  with open("validation_files/valRK4_time_Ni{}.txt".format(N)) as s:
    lines = s.read().splitlines()
    #for line in lines:
    temp=lines[0].split()
    for numb in temp:
      templist.append(float(numb))

  time_RK4[N] = np.array(templist)

  plt.figure(figsize=(h_size, v_size))
  plt.plot(time_RK4[N]*10**(6),solution_exact_RK4[N])
  plt.plot(time_RK4[N]*10**(6),solution_linear_RK4[N])
  plt.plot(time_RK4[N]*10**(6),solution_nonlinear_RK4[N])
  plt.title( "$N_i$ = {}     $n1_0$/$n_0$ = {}".format(N,n10_n0) )
  plt.ylabel("$n_1$/$n_0$")
  plt.xlabel("time [$\mu s$]")
  plt.legend(["exact","linear","nonlinear"])
  plt.tight_layout() 
  filename = "validation_files/solution_Ni{}_{}".format(N,n10_n0)
  filename = filename.replace(".","_") + ".png"
  plt.savefig(filename)

for N in Ni_v:
  max_index = np.argmax( np.abs (solution_linear_RK4[N] - solution_exact_RK4[N])  )
  err_RK4.append( np.abs( ( solution_linear_RK4[N][max_index] - solution_exact_RK4[N][max_index]  ) \
                   / solution_exact_RK4[N][max_index] ) )

Ni_v = np.array(Ni_v)
err_RK4 = np.array(err_RK4)

plt.figure(figsize=(h_size, v_size))
plt.loglog(Ni_v,err_RK4,'*')
plt.title("$n1_0$/$n_0$ = {}".format(n10_n0))
plt.xlabel("$N_i$")
plt.ylabel("error")
plt.tight_layout() 
plt.savefig("validation_files/Convergence_RK4.png") 

plt.show()

