import matplotlib.font_manager
matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
formatter = ticker.ScalarFormatter(useMathText=True) #scientific notation
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 

font = {'family' : 'Helvetica',
        'size'   : 22}

		
plt.rc('font', **font)
plt.rcParams['axes.unicode_minus']=False


h_size = 9
v_size = 6
linesize = 3.0
linesize_large = 4.0

import numpy as np
import os
import sys

# print(len(sys.argv))

if len(sys.argv)==1:
  solution_file = "solution.txt"
  isexperiment = 0
elif len(sys.argv)==2:
  solution_file=sys.argv[1]
  isexperiment=0
elif len(sys.argv)==3:
  solution_file=sys.argv[1]
  data_file=sys.argv[2]
  isexperiment=1
else:
  print("Wrong number of inputs. Usage: python3 plotsolution.py <solution_file> <data_file> (optional)")
  sys.exit()



with open(solution_file) as s:
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
        # mode independent values
        plasma_potential = np.zeros(N_voltages, dtype = np.double)
        ion_velocity0 = np.zeros(N_voltages, dtype = np.double)
        electron_diffusion = np.zeros(N_voltages, dtype = np.double)
        electron_temperature = np.zeros(N_voltages, dtype = np.double)
        Jxi0 = np.zeros(N_voltages, dtype = np.double)
        # decreasing modes values
        m_decr = np.zeros(N_voltages, dtype = int)	
        frequency_decr = np.zeros(N_voltages, dtype = np.double)
        dphi_decr = np.zeros(N_voltages, dtype = np.double)
        vi1_decr = np.zeros(N_voltages, dtype = np.double)
        ve1_decr = np.zeros(N_voltages, dtype = np.double)
        dJi1_decr = np.zeros(N_voltages, dtype = np.double)
        dJe1_decr = np.zeros(N_voltages, dtype = np.double)
        dIix_decr = np.zeros(N_voltages, dtype = np.double)
        dIex_decr = np.zeros(N_voltages, dtype = np.double)
        Ix_decr = np.zeros(N_voltages, dtype = np.double)
        Lbmin_decr = np.zeros(N_voltages, dtype = np.double)
        A_decr = np.zeros(N_voltages, dtype = np.double)
        B_decr = np.zeros(N_voltages, dtype = np.double)
        C_decr = np.zeros(N_voltages, dtype = np.double)
        alpha_decr = np.zeros(N_voltages, dtype = np.double)
        beta_decr = np.zeros(N_voltages, dtype = np.double)
        gamma_decr = np.zeros(N_voltages, dtype = np.double)
        # increasing modes values
        m_incr = np.zeros(N_voltages, dtype = int)	
        frequency_incr = np.zeros(N_voltages, dtype = np.double)
        dphi_incr = np.zeros(N_voltages, dtype = np.double)
        vi1_incr = np.zeros(N_voltages, dtype = np.double)
        ve1_incr = np.zeros(N_voltages, dtype = np.double)
        dJi1_incr = np.zeros(N_voltages, dtype = np.double)
        dJe1_incr = np.zeros(N_voltages, dtype = np.double)
        dIix_incr = np.zeros(N_voltages, dtype = np.double)
        dIex_incr = np.zeros(N_voltages, dtype = np.double)
        Ix_incr = np.zeros(N_voltages, dtype = np.double)
        Lbmin_incr = np.zeros(N_voltages, dtype = np.double)
        A_incr = np.zeros(N_voltages, dtype = np.double)
        B_incr = np.zeros(N_voltages, dtype = np.double)
        C_incr = np.zeros(N_voltages, dtype = np.double)
        alpha_incr = np.zeros(N_voltages, dtype = np.double)
        beta_incr = np.zeros(N_voltages, dtype = np.double)
        gamma_incr = np.zeros(N_voltages, dtype = np.double)
    
        i = 0
      else:
        try:
          # mode independent values
          temp_value = np.double(temp[0])
          plasma_potential[i] = temp_value
          ion_velocity0[i] = np.double(temp[1])
          electron_diffusion[i] = np.double(temp[2])
          electron_temperature[i] = np.double(temp[3])
          Jxi0[i] = np.double(temp[4])
          # decreasing modes values
          m_decr[i] = np.double(temp[5])
          frequency_decr[i] = np.double(temp[6])
          dphi_decr[i] = np.double(temp[7])
          vi1_decr[i] = np.double(temp[8])
          ve1_decr[i] = np.double(temp[9])
          dJi1_decr[i] = np.double(temp[10])
          dJe1_decr[i] = np.double(temp[11])
          dIix_decr[i] = np.double(temp[12])
          dIex_decr[i] = np.double(temp[13])
          Ix_decr[i] = np.double(temp[14])
          Lbmin_decr[i] = np.double(temp[15])
          A_decr[i] = np.double(temp[16])
          B_decr[i] = np.double(temp[17])
          C_decr[i] = np.double(temp[18])
          alpha_decr[i] = np.double(temp[19])
          beta_decr[i] = np.double(temp[20])
          gamma_decr[i] = np.double(temp[21])
          # increasing modes values
          m_incr[i] = np.double(temp[22])
          frequency_incr[i] = np.double(temp[23])
          dphi_incr[i] = np.double(temp[24])
          vi1_incr[i] = np.double(temp[25])
          ve1_incr[i] = np.double(temp[26])
          dJi1_incr[i] = np.double(temp[27])
          dJe1_incr[i] = np.double(temp[28])
          dIix_incr[i] = np.double(temp[29])
          dIex_incr[i] = np.double(temp[30])
          Ix_incr[i] = np.double(temp[31])
          Lbmin_incr[i] = np.double(temp[32])
          A_incr[i] = np.double(temp[33])
          B_incr[i] = np.double(temp[34])
          C_incr[i] = np.double(temp[35])
          alpha_incr[i] = np.double(temp[36])
          beta_incr[i] = np.double(temp[37])
          gamma_incr[i] = np.double(temp[38])
          i = i + 1
        except(ValueError):
          pass
    except(IndexError):
      pass
  
if isexperiment: 
  with open(data_file) as s:
    lines = s.read().splitlines()
    gap_voltage_exp=np.zeros(len(lines)-2)
    gap_current_exp=np.zeros(len(lines)-2)
    mode_number_exp=np.zeros(len(lines)-2)
    frequency_exp=np.zeros(len(lines)-2)
    for i in range(1,len(gap_voltage_exp)+1):
      line = lines[i]
      temp=line.split()
      try:
        gap_voltage_exp[i-1]=np.double(temp[0])
        gap_current_exp[i-1]=np.double(temp[1])
        mode_number_exp[i-1]=np.double(temp[2])
        frequency_exp[i-1]=np.double(temp[3])
      except(ValueError):
        pass
      except(IndexError):
        pass
    index_min = np.argmin(gap_voltage_exp)
    gap_voltage_exp_decr=gap_voltage_exp[:index_min+1]
    gap_current_exp_decr=gap_current_exp[:index_min+1]
    mode_number_exp_decr=mode_number_exp[:index_min+1]
    frequency_exp_decr=frequency_exp[:index_min+1]
    gap_voltage_exp_incr=gap_voltage_exp[index_min+1:]
    gap_current_exp_incr=gap_current_exp[index_min+1:]
    mode_number_exp_incr=mode_number_exp[index_min+1:]
    frequency_exp_incr=frequency_exp[index_min+1:]				
    voltage_shift = np.double(lines[-1].split()[0])

	
#unique_m = set(list(set(m_decr)) + list(set(m_incr)))
unique_m = (2, 3, 4, 5, 6, 7, 8)
legend_list=[]
colors = {np.nan: 'k', 2: (0.6350, 0.0780, 0.1840),3: (0, 0.4470, 0.7410),4: (0.8500, 0.3250, 0.0980), 5: (0.9290, 0.6940, 0.1250), 6: (0.4940, 0.1840, 0.5560), 1: (0.4660, 0.6740, 0.1880), 7 : (0.3010, 0.7450, 0.9330), 8 : (1 , 0, 0)}

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, solution_file[:-4] + '_figures/')
if not os.path.isdir(results_dir):
  os.makedirs(results_dir)
	
plt.figure(1,figsize=(h_size, v_size))
for mode in unique_m:
  if mode in m_decr:
    plt.plot(plasma_potential[m_decr == mode],frequency_decr[m_decr == mode]/1000,color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)]                                       
  if mode in m_incr:                                                                                  
    plt.plot(plasma_potential[m_incr == mode],frequency_incr[m_incr == mode]/1000,'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]                                       

if isexperiment:
  for mode in set(mode_number_exp_decr[~np.isnan(mode_number_exp_decr)]):
    plt.plot(abs(gap_voltage_exp_decr[mode_number_exp_decr == mode])-voltage_shift,frequency_exp_decr[mode_number_exp_decr == mode]/1000,'v',color = colors[mode])	
    legend_list = legend_list + ["m$_{decr}$ = " + str(int(mode))]
  for mode in set(mode_number_exp_incr[~np.isnan(mode_number_exp_incr)]):
    plt.plot(abs(gap_voltage_exp_incr[mode_number_exp_incr == mode])-voltage_shift,frequency_exp_incr[mode_number_exp_incr == mode]/1000,'^',color = colors[mode])	
    legend_list = legend_list + ["m$_{incr}$ = " + str(int(mode))]
#Axes.ticklabel_format(self, *, axis='both', style='', scilimits=None, useOffset=None, useLocale=None, useMathText=None)
	
#ax1.set_aspect('equal')	
#plt.axes().set_aspect('equal')
plt.legend(legend_list,loc="upper right" , fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Frequency [kHz]")
plt.tight_layout()
#plt.set_aspect('equal')
#plt.show(block=True)
plt.savefig(results_dir+"f_vs_EL.png")


plt.figure(2,figsize=(h_size, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_decr:
    plt.plot(plasma_potential[m_decr == mode],Ix_decr[m_decr == mode]*1000,color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
  if mode in m_incr:
    plt.plot(plasma_potential[m_incr == mode],Ix_incr[m_incr == mode]*1000,'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]

if isexperiment:
  for mode in set(mode_number_exp_decr[~np.isnan(mode_number_exp_decr)]):
    plt.plot(abs(gap_voltage_exp_decr[mode_number_exp_decr == mode])-voltage_shift,gap_current_exp_decr[mode_number_exp_decr == mode],'v',color = colors[mode])	
    legend_list = legend_list + ["m$_{decr}$ = " + str(int(mode))]
  for mode in set(mode_number_exp_incr[~np.isnan(mode_number_exp_incr)]):
    plt.plot(abs(gap_voltage_exp_incr[mode_number_exp_incr == mode])-voltage_shift,gap_current_exp_incr[mode_number_exp_incr == mode],'^',color = colors[mode])	
    legend_list = legend_list + ["m$_{incr}$ = " + str(int(mode))]


plt.legend(legend_list, fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Current [mA]")
plt.tight_layout()
#plt.show(block=True)
plt.savefig(results_dir+"I_vs_EL.png")
plt.savefig(results_dir+"I_vs_EL.eps")

plt.figure(3,figsize=(h_size, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_decr:
    plt.plot(plasma_potential[m_decr == mode],dphi_decr[m_decr == mode],color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
  if mode in m_incr:
    plt.plot(plasma_potential[m_incr == mode],dphi_incr[m_incr == mode],'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]

plt.legend(legend_list , fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("$d\phi$ [V]")
plt.tight_layout()
#plt.show(block=True)
plt.savefig(results_dir+"ne_over_n0_vs_EL.png")
plt.savefig(results_dir+"ne_over_n0_vs_EL.eps")

plt.figure(4,figsize=(h_size, v_size))
legend_list=[]

for mode in unique_m:
  if mode in m_decr:
    plt.plot(plasma_potential[m_decr == mode],vi1_decr[m_decr == mode],color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
  if mode in m_incr:
    plt.plot(plasma_potential[m_incr == mode],vi1_incr[m_incr == mode],'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]

plt.legend(legend_list, fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("v$_{i1}$ [m/s]")
#plt.show(block=True)
plt.tight_layout()
plt.savefig(results_dir+"vi1_vs_EL.png")
plt.savefig(results_dir+"vi1_vs_EL.eps")    
	
	
plt.figure(5,figsize=(h_size*1.05, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_decr:
    plt.plot(plasma_potential[m_decr == mode],frequency_decr[m_decr == mode]*2*np.pi,color=colors[mode], linewidth = linesize_large)
    legend_list = legend_list + ["m = " + str(mode)] 
    plt.plot(plasma_potential[m_decr == mode],A_decr[m_decr == mode],'--',color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["A"] 
    plt.plot(plasma_potential[m_decr == mode],B_decr[m_decr == mode],'-.',color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["B"] 
    plt.plot(plasma_potential[m_decr == mode],C_decr[m_decr == mode],':',color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["C"]

		
		
art = plt.legend(legend_list, fancybox = False, edgecolor = 'k' , loc='right' , bbox_to_anchor=(1.425, 0.5) , fontsize = 22, labelspacing=0.188)
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("$\omega_R$ [rad/s]")

#plt.gca().yaxis.set_minor_formatter(ticker.NullFormatter())
#plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter())
#plt.gca().yaxis.get_major_formatter().set_scientific(False)
#plt.gca().yaxis.get_major_formatter().set_useOffset(False)

plt.tight_layout()
#plt.show(block=True)
plt.savefig(results_dir+"wR_vs_EL_decr_comparison.png", additional_artists = art, bbox_inches="tight")
plt.savefig(results_dir+"wR_vs_EL_decr_comparison.eps", additional_artists = art, bbox_inches="tight")


plt.figure(6,figsize=(h_size, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_incr:
    plt.plot(plasma_potential[m_incr == mode],frequency_incr[m_incr == mode]*2*np.pi,color=colors[mode], linewidth = linesize_large)
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)] 
    plt.plot(plasma_potential[m_incr == mode],A_incr[m_incr == mode],'--',color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["A"] 
    plt.plot(plasma_potential[m_incr == mode],B_incr[m_incr == mode],'-.',color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["B"] 
    plt.plot(plasma_potential[m_incr == mode],C_incr[m_incr == mode],':',color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["C"] 



art = plt.legend(legend_list, loc='right', fancybox = False, edgecolor = 'k' , bbox_to_anchor=(1.4, 0.5), fontsize = 22, labelspacing=0.07)
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("$\omega_R$ [rad/s]")
plt.tight_layout()
#plt.show(block=True)
plt.savefig(results_dir+"wR_vs_EL_incr_comparison.png", additional_artists = art, bbox_inches="tight")
plt.savefig(results_dir+"wR_vs_EL_incr_comparison.eps", additional_artists = art, bbox_inches="tight")


#                 W                       #


with open("w_solution.txt") as s:
  lines = s.read().splitlines()
  inputs = ("ion_mass_u", "Te", "B0", "R0", "LB", "Ln", "kz", "kx", "N_voltages")
  
  ky = []
  wR = []
  A = []
  B = []
  C = []
  wI = []
  alpha = []
  beta = []
  gamma = []
  delta = []
  
  for line in lines:
    # first I read the input parameters
    # I need to know N_voltages to initialize
    # the numpy arrays.
    temp=line.split()
    
    ky.append(np.double(temp[0])/R0) 
    wR.append(np.double(temp[1])) 
    A.append(np.double(temp[2])) 
    B.append(np.double(temp[3])) 
    C.append(np.double(temp[4])) 
    wI.append(np.double(temp[5]))
    alpha.append(np.double(temp[6])) 
    beta.append(np.double(temp[7])) 
    gamma.append(np.double(temp[8]))

  ky = np.array(ky, dtype=np.double) 
  wR = np.array(wR, dtype=np.double) 
  A = np.array(A, dtype=np.double) 
  B = np.array(B, dtype=np.double) 
  C = np.array(C, dtype=np.double) 
  wI = np.array(wI, dtype=np.double) 
  alpha = np.array(alpha, dtype=np.double) 
  beta = np.array(beta, dtype=np.double) 
  gamma = np.array(gamma, dtype=np.double) 
  
  plt.figure(7,figsize=(h_size, v_size))
  plt.plot(ky*R0,wR, linewidth = linesize_large)	
  plt.plot(ky*R0,A,'--', linewidth = linesize)	
  plt.plot(ky*R0,B,'-.', linewidth = linesize)	
  plt.plot(ky*R0,C,':', linewidth = linesize)	
  plt.legend(["$\omega_R$ tot", "A", "B", "C"], fancybox = False, edgecolor = 'k' , loc = 'upper right', ncol=2)
  #pylab.legend(loc=9, bbox_to_anchor=(0.5, -0.1))
  plt.xlabel("m = k$_y$R$_0$")
  plt.ylabel("$\omega{_R}$ [rad/s]")
  plt.tight_layout()
  
  plt.savefig(results_dir+"wR_vs_ky.png")
  plt.savefig(results_dir+"wR_vs_ky.eps")
  
  plt.figure(8,figsize=(h_size, v_size))
  

  plt.plot(ky*R0,wI*wI, linewidth = linesize)
  plt.plot(ky*R0,alpha, linewidth = linesize)
  plt.plot(ky*R0,beta, linewidth = linesize)
  plt.plot(ky*R0,gamma, linewidth = linesize)
  plt.legend(["$\omega{_I}^2$", "$α$", "$β$","$γ$"], loc = 'upper right' , fancybox = False, edgecolor = 'k' )
  plt.xlabel("m = k$_y$R$_0$")
  plt.ylabel("[rad/s]$^2$")
  plt.savefig(results_dir+"wI_vs_ky.png")
  plt.savefig(results_dir+"wI_vs_ky.eps")
  
  plt.figure(9,figsize=(h_size, v_size))
  plt.plot(ky*R0,[alpha[i]+beta[i]+gamma[i] for i in range(len(alpha))],'r', linewidth = linesize_large)
  plt.plot(ky*R0,[alpha[i]-beta[i]-gamma[i] for i in range(len(alpha))],'b', linewidth = linesize_large)
  plt.legend(["C$^2$", "$\omega{_I}^2$"], loc = 'upper right', fancybox = False, edgecolor = 'k' )
  plt.xlabel("m = k$_y$R$_0$")
  plt.ylabel("[rad/s]$^2$")
  plt.tight_layout()
  
  plt.savefig(results_dir+"sqrt_vs_ky.png")
  plt.savefig(results_dir+"sqrt_vs_ky.eps")
  
  plt.figure(10,figsize=(h_size, v_size))
  #plt.plot(ky,[alpha[i]-beta[i]-gamma[i] for i in range(len(alpha))],'b', linewidth = linesize_large)
  #plt.plot(ky,[alpha[i]+beta[i]+gamma[i] for i in range(len(alpha))],'r', linewidth = linesize)
  plt.plot(ky*R0,[beta[i]+gamma[i] for i in range(len(alpha))], linewidth = linesize)
  plt.plot(ky*R0,alpha,'--', linewidth = linesize)	
  plt.plot(ky*R0,beta,'-.', linewidth = linesize)	
  plt.plot(ky*R0,gamma,':', linewidth = linesize)	
  leg = plt.legend(["$β + γ$", "$α$", "$β$","$γ$"], loc = 'upper right', fontsize = 20 , ncol=2, fancybox = False, edgecolor = 'k' , handleheight=2.0, labelspacing=0.04)
  plt.xlabel("m = k$_y$R$_0$")
  plt.ylabel("[rad/s]$^2$")
  plt.tight_layout()
  
  plt.savefig(results_dir+"sqrt_vs_ky.png")
  plt.savefig(results_dir+"sqrt_vs_ky.eps")
  
  plt.show()

