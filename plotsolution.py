import matplotlib.font_manager
matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.animation as animation
import scipy.integrate as integrate

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

# subset of modes
unique_m = (3, 4)
do_decr = True
do_incr = True
do_movies = False

import numpy as np
import os
import sys

# print(len(sys.argv))

if len(sys.argv)==1:
  solution_file = "linear_model_solutions/solution.txt"
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
        # mode independent values
        plasma_potential = np.zeros(N_voltages, dtype = np.double)
        ne0 = np.zeros(N_voltages, dtype = np.double)
        ion_velocity0 = np.zeros(N_voltages, dtype = np.double)
        electron_diffusion = np.zeros(N_voltages, dtype = np.double)
        electron_temperature = np.zeros(N_voltages, dtype = np.double)
        Jxi0 = np.zeros(N_voltages, dtype = np.double)
        Ixi0 = np.zeros(N_voltages, dtype = np.double)
        # decreasing modes values
        m_decr = np.zeros(N_voltages, dtype = int)	
        frequency_decr = np.zeros(N_voltages, dtype = np.double)
        dphi_decr = np.zeros(N_voltages, dtype = np.double)
        vi1_decr = np.zeros(N_voltages, dtype = np.double)
        ve1_EXB_decr = np.zeros(N_voltages, dtype = np.double)
        ve1_D_decr = np.zeros(N_voltages, dtype = np.double)
        
        dphi_phase_decr = np.zeros(N_voltages, dtype = np.double)
        vi1_phase_decr = np.zeros(N_voltages, dtype = np.double)
        ve1_EXB_phase_decr = np.zeros(N_voltages, dtype = np.double)
        ve1_D_phase_decr = np.zeros(N_voltages, dtype = np.double)
        
        dJi1_ave_decr = np.zeros(N_voltages, dtype = np.double)
        dJe1_EXB_ave_decr = np.zeros(N_voltages, dtype = np.double)
        dJe1_D_ave_decr = np.zeros(N_voltages, dtype = np.double)       
        dIix_ave_decr = np.zeros(N_voltages, dtype = np.double)
        dIex_EXB_ave_decr = np.zeros(N_voltages, dtype = np.double)
        dIex_D_ave_decr = np.zeros(N_voltages, dtype = np.double)
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
        ve1_EXB_incr = np.zeros(N_voltages, dtype = np.double)
        ve1_D_incr = np.zeros(N_voltages, dtype = np.double)
        
        dphi_phase_incr = np.zeros(N_voltages, dtype = np.double)
        vi1_phase_incr = np.zeros(N_voltages, dtype = np.double)
        ve1_EXB_phase_incr = np.zeros(N_voltages, dtype = np.double)
        ve1_D_phase_incr = np.zeros(N_voltages, dtype = np.double)
        
        dJi1_ave_incr = np.zeros(N_voltages, dtype = np.double)
        dJe1_EXB_ave_incr = np.zeros(N_voltages, dtype = np.double)
        dJe1_D_ave_incr = np.zeros(N_voltages, dtype = np.double)
        dIix_ave_incr = np.zeros(N_voltages, dtype = np.double)
        dIex_EXB_ave_incr = np.zeros(N_voltages, dtype = np.double)
        dIex_D_ave_incr = np.zeros(N_voltages, dtype = np.double)
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
          input_number=0;
          temp_value = np.double(temp[input_number])
          plasma_potential[i] = temp_value
          input_number+=1
          ne0[i] = np.double(temp[input_number])
          input_number+=1          
          ion_velocity0[i] = np.double(temp[input_number])
          input_number+=1
          electron_diffusion[i] = np.double(temp[input_number])
          input_number+=1
          electron_temperature[i] = np.double(temp[input_number])
          input_number+=1
          Jxi0[i] = np.double(temp[input_number])
          input_number+=1
          Ixi0[i] = np.double(temp[input_number])
          input_number+=1
          # decreasing modes values
          m_decr[i] = np.double(temp[input_number])
          input_number+=1
          frequency_decr[i] = np.double(temp[input_number])
          input_number+=1
          dphi_decr[i] = np.double(temp[input_number])
          input_number+=1
          vi1_decr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_EXB_decr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_D_decr[i] = np.double(temp[input_number])
          input_number+=1
          
          dphi_phase_decr[i] = np.double(temp[input_number])
          input_number+=1
          vi1_phase_decr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_EXB_phase_decr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_D_phase_decr[i] = np.double(temp[input_number])
          input_number+=1
          
          dJi1_ave_decr[i] = np.double(temp[input_number])
          input_number+=1
          dJe1_EXB_ave_decr[i] = np.double(temp[input_number])
          input_number+=1
          dJe1_D_ave_decr[i] = np.double(temp[input_number])
          input_number+=1
          dIix_ave_decr[i] = np.double(temp[input_number])
          input_number+=1
          dIex_EXB_ave_decr[i] = np.double(temp[input_number])
          input_number+=1
          dIex_D_ave_decr[i] = np.double(temp[input_number])
          input_number+=1       
          Ix_decr[i] = np.double(temp[input_number])
          input_number+=1
          Lbmin_decr[i] = np.double(temp[input_number])
          input_number+=1
          A_decr[i] = np.double(temp[input_number])
          input_number+=1
          B_decr[i] = np.double(temp[input_number])
          input_number+=1
          C_decr[i] = np.double(temp[input_number])
          input_number+=1
          alpha_decr[i] = np.double(temp[input_number])
          input_number+=1
          beta_decr[i] = np.double(temp[input_number])
          input_number+=1
          gamma_decr[i] = np.double(temp[input_number])
          input_number+=1
          # increasing modes values
          m_incr[i] = np.double(temp[input_number])
          input_number+=1
          frequency_incr[i] = np.double(temp[input_number])
          input_number+=1
          dphi_incr[i] = np.double(temp[input_number])
          input_number+=1
          vi1_incr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_EXB_incr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_D_incr[i] = np.double(temp[input_number])
          input_number+=1
          
          dphi_phase_incr[i] = np.double(temp[input_number])
          input_number+=1
          vi1_phase_incr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_EXB_phase_incr[i] = np.double(temp[input_number])
          input_number+=1
          ve1_D_phase_incr[i] = np.double(temp[input_number])
          input_number+=1
 
          dJi1_ave_incr[i] = np.double(temp[input_number])
          input_number+=1
          dJe1_EXB_ave_incr[i] = np.double(temp[input_number])
          input_number+=1
          dJe1_D_ave_incr[i] = np.double(temp[input_number])
          input_number+=1
          dIix_ave_incr[i] = np.double(temp[input_number])
          input_number+=1
          dIex_EXB_ave_incr[i] = np.double(temp[input_number])
          input_number+=1
          dIex_D_ave_incr[i] = np.double(temp[input_number])
          input_number+=1
          Ix_incr[i] = np.double(temp[input_number])
          input_number+=1
          Lbmin_incr[i] = np.double(temp[input_number])
          input_number+=1
          A_incr[i] = np.double(temp[input_number])
          input_number+=1
          B_incr[i] = np.double(temp[input_number])
          input_number+=1
          C_incr[i] = np.double(temp[input_number])
          input_number+=1
          alpha_incr[i] = np.double(temp[input_number])
          input_number+=1
          beta_incr[i] = np.double(temp[input_number])
          input_number+=1
          gamma_incr[i] = np.double(temp[input_number])
          input_number+=1
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

	
legend_list=[]
colors = {np.nan: 'k', 2: (0.6350, 0.0780, 0.1840),3: (0, 0.4470, 0.7410),4: (0.8500, 0.3250, 0.0980), 5: (0.9290, 0.6940, 0.1250), 6: (0.4940, 0.1840, 0.5560), 1: (0.4660, 0.6740, 0.1880), 7 : (0.3010, 0.7450, 0.9330), 8 : (1 , 0, 0)}

script_dir = os.path.dirname(__file__)
figures_dir = 'solution_figures/'
movies_dir = 'solution_movies/'

if not os.path.isdir(figures_dir):
  os.makedirs(figures_dir)
 
if not os.path.isdir(movies_dir):
  os.makedirs(movies_dir) 
	
plt.figure(1,figsize=(h_size, v_size))
for mode in unique_m:
  if mode in m_decr and do_decr:
    plt.plot(plasma_potential[m_decr == mode],frequency_decr[m_decr == mode]/1000,color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)]                                       
  if mode in m_incr and do_incr:                                                                                  
    plt.plot(plasma_potential[m_incr == mode],frequency_incr[m_incr == mode]/1000,'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]                                       

if isexperiment:
  if do_decr:
    for mode in set(mode_number_exp_decr[~np.isnan(mode_number_exp_decr)]):
      plt.plot(abs(gap_voltage_exp_decr[mode_number_exp_decr == mode])-voltage_shift,frequency_exp_decr[mode_number_exp_decr == mode]/1000,'v',markersize=10,color = colors[mode])	
      legend_list = legend_list + ["m$_{decr}$ = " + str(int(mode))]
  if do_incr:
    for mode in set(mode_number_exp_incr[~np.isnan(mode_number_exp_incr)]):
      plt.plot(abs(gap_voltage_exp_incr[mode_number_exp_incr == mode])-voltage_shift,frequency_exp_incr[mode_number_exp_incr == mode]/1000,'^',markersize=10,color = colors[mode])	
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
plt.savefig(figures_dir+"f_vs_EL.png")


plt.figure(2,figsize=(h_size, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_decr and do_decr:
    plt.plot(plasma_potential[m_decr == mode],Ix_decr[m_decr == mode]*1000,color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
  if mode in m_incr and do_incr:
    plt.plot(plasma_potential[m_incr == mode],Ix_incr[m_incr == mode]*1000,'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]

if isexperiment:
  if do_decr:
    for mode in set(mode_number_exp_decr[~np.isnan(mode_number_exp_decr)]):
      plt.plot(abs(gap_voltage_exp_decr[mode_number_exp_decr == mode])-voltage_shift,gap_current_exp_decr[mode_number_exp_decr == mode],'v',markersize=10,color = colors[mode])	
      legend_list = legend_list + ["m$_{decr}$ = " + str(int(mode))]
  if do_incr:
    for mode in set(mode_number_exp_incr[~np.isnan(mode_number_exp_incr)]):
      plt.plot(abs(gap_voltage_exp_incr[mode_number_exp_incr == mode])-voltage_shift,gap_current_exp_incr[mode_number_exp_incr == mode],'^',markersize=10,color = colors[mode])	
      legend_list = legend_list + ["m$_{incr}$ = " + str(int(mode))]


plt.legend(legend_list, fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Current [mA]")
plt.tight_layout()
#plt.show(block=True)
plt.savefig(figures_dir+"I_vs_EL.png")
plt.savefig(figures_dir+"I_vs_EL.eps")



plt.figure(20,figsize=(h_size, v_size))
legend_list=[]
plt.plot(plasma_potential,Ixi0*1000,"k", linewidth = linesize)
legend_list = legend_list + ["$I_{ix_0}$"]     
for mode in unique_m:
  if mode in m_decr and do_decr:
    plt.plot(plasma_potential[m_decr == mode],dIex_EXB_ave_decr[m_decr == mode]*1000,'--',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$I_{ex_{ExB}}$"]     
    plt.plot(plasma_potential[m_decr == mode],dIex_D_ave_decr[m_decr == mode]*1000,'-.',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$I_{ex_{D}}$"] 
    plt.plot(plasma_potential[m_decr == mode],dIix_ave_decr[m_decr == mode]*1000,':',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$I_{ix}$"]         
    plt.plot(plasma_potential[m_decr == mode],Ix_decr[m_decr == mode]*1000,color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["$I_{x_{tot}}$"] 
  if mode in m_incr and do_incr:
    plt.plot(plasma_potential[m_decr == mode],dIex_EXB_ave_incr[m_decr == mode]*1000,'--',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$I_{ex_{ExB}}$"]     
    plt.plot(plasma_potential[m_decr == mode],dIex_D_ave_incr[m_decr == mode]*1000,'-.',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$I_{ex_{D}}$"]  
    plt.plot(plasma_potential[m_decr == mode],dIix_ave_incr[m_decr == mode]*1000,':',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$I_{ix}$"]   
    plt.plot(plasma_potential[m_incr == mode],Ix_incr[m_incr == mode]*1000,color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["$I_{x_{tot}}$"]

if isexperiment:
  if do_decr:
    for mode in set(mode_number_exp_decr[~np.isnan(mode_number_exp_decr)]):
      plt.plot(abs(gap_voltage_exp_decr[mode_number_exp_decr == mode])-voltage_shift,gap_current_exp_decr[mode_number_exp_decr == mode],'v',markersize=10,color = colors[mode])	
      legend_list = legend_list + ["m$_{decr}$ = " + str(int(mode))]
  if do_incr:
    for mode in set(mode_number_exp_incr[~np.isnan(mode_number_exp_incr)]):
      plt.plot(abs(gap_voltage_exp_incr[mode_number_exp_incr == mode])-voltage_shift,gap_current_exp_incr[mode_number_exp_incr == mode],'^',markersize=10,color = colors[mode])	
      legend_list = legend_list + ["m$_{incr}$ = " + str(int(mode))]


plt.legend(legend_list, fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Current [mA]")
plt.tight_layout()
#plt.show(block=True)
plt.savefig(figures_dir+"I_vs_EL.png")
plt.savefig(figures_dir+"I_vs_EL.eps")


# 
plt.figure(21,figsize=(h_size, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_decr and do_decr:
    plt.plot(plasma_potential[m_decr == mode],dphi_phase_decr[m_decr == mode],'--',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$\phi$"]     
    plt.plot(plasma_potential[m_decr == mode],vi1_phase_decr[m_decr == mode],'-.',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$v_i$"] 
    plt.plot(plasma_potential[m_decr == mode],ve1_EXB_phase_decr[m_decr == mode],':',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$v_e{_ExB}$"]         
    plt.plot(plasma_potential[m_decr == mode],ve1_D_phase_decr[m_decr == mode],color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["$v_e{_D}$"] 
  if mode in m_incr and do_incr:
    plt.plot(plasma_potential[m_decr == mode],dphi_phase_decr[m_decr == mode],'--',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$\phi$"]     
    plt.plot(plasma_potential[m_decr == mode],vi1_phase_decr[m_decr == mode],'-.',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$v_i$"] 
    plt.plot(plasma_potential[m_decr == mode],ve1_EXB_phase_decr[m_decr == mode],':',color=colors[mode], linewidth = linesize)    
    legend_list = legend_list + ["$v_e{_ExB}$"]         
    plt.plot(plasma_potential[m_decr == mode],ve1_D_phase_decr[m_decr == mode],color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["$v_e{_D}$"] 

plt.legend(legend_list, fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Phase [rad]")
plt.tight_layout()
#plt.show(block=True)
figure_name = "phase_shift"
plt.savefig(figures_dir+figure_name+".png")
plt.savefig(figures_dir+figure_name+".eps")


plt.figure(3,figsize=(h_size, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_decr and do_decr:
    plt.plot(plasma_potential[m_decr == mode],dphi_decr[m_decr == mode],color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
  if mode in m_incr and do_incr:
    plt.plot(plasma_potential[m_incr == mode],dphi_incr[m_incr == mode],'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]

plt.legend(legend_list , fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("$d\phi$ [V]")
plt.tight_layout()
#plt.show(block=True)
plt.savefig(figures_dir+"ne_over_n0_vs_EL.png")
plt.savefig(figures_dir+"ne_over_n0_vs_EL.eps")

plt.figure(4,figsize=(h_size, v_size))
legend_list=[]

for mode in unique_m:
  if mode in m_decr and do_decr:
    plt.plot(plasma_potential[m_decr == mode],vi1_decr[m_decr == mode],color=colors[mode], linewidth = linesize)
    legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
  if mode in m_incr and do_incr:
    plt.plot(plasma_potential[m_incr == mode],vi1_incr[m_incr == mode],'--',color = colors[mode], linewidth = linesize)	
    legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]

plt.legend(legend_list, fancybox = False, edgecolor = 'k' )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("v$_{i1}$ [m/s]")
#plt.show(block=True)
plt.tight_layout()
plt.savefig(figures_dir+"vi1_vs_EL.png")
plt.savefig(figures_dir+"vi1_vs_EL.eps")    
	
	
plt.figure(5,figsize=(h_size*1.05, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_decr and do_decr:
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
plt.savefig(figures_dir+"wR_vs_EL_decr_comparison.png", additional_artists = art, bbox_inches="tight")
plt.savefig(figures_dir+"wR_vs_EL_decr_comparison.eps", additional_artists = art, bbox_inches="tight")


plt.figure(6,figsize=(h_size, v_size))
legend_list=[]
for mode in unique_m:
  if mode in m_incr and do_incr:
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
plt.savefig(figures_dir+"wR_vs_EL_incr_comparison.png", additional_artists = art, bbox_inches="tight")
plt.savefig(figures_dir+"wR_vs_EL_incr_comparison.eps", additional_artists = art, bbox_inches="tight")




# MOVIES

def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



# 0.568743859261
e = 1.60217657e-19

Lz = 0.001

def animate_currents(ind):
  # ind = find_nearest_index(ion_velocity0,-500)
  m_decr_ = m_decr[ind]
  ky_ = m_decr_/R0
  f_  = frequency_decr[ind]
  lambda_ = 2*np.pi/ky_
  y_ = np.linspace(0,m_decr_*lambda_,1000)
  ne0_ = ne0[ind]*(1.0+np.sin(ky_*y_))
  De_ = electron_diffusion[ind]

  dphi_phase_ = dphi_phase_decr[ind]
  vi1_phase_ = vi1_phase_decr[ind]
  ve1_EXB_phase_ = ve1_EXB_phase_decr[ind] 
  ve1_D_phase_ = ve1_D_phase_decr[ind]

  dphi_ = plasma_potential[ind] + dphi_decr[ind] * np.sin(ky_*y_ + dphi_phase_)
  vi0_ = ion_velocity0[ind]
  vi1_ = vi1_decr[ind]*np.sin(ky_*y_ + vi1_phase_)
  ve1_EXB_ = ve1_EXB_decr[ind]*np.sin(ky_*y_ + ve1_EXB_phase_)
  ve1_D_ = ve1_D_decr[ind]*np.sin(ky_*y_ + ve1_D_phase_)
  
  ve_De_ = De_/Lz  + kz*De_*np.sin(ky_*y_ + np.pi/2)
  

  Je_De_ = -e*ne0_*ve_De_
  Ji0_ = e*ne0_*vi0_
  Ji1_ =  e*ne0_*vi1_
  JeExB_ = -e*ne0_*ve1_EXB_
  JeD_ = -e*ne0_*ve1_D_
  Jtot_ = Ji0_ + Ji1_ + JeExB_ + JeD_ + Je_De_

  # plt.figure(40,figsize=(h_size, v_size))
  # plt.plot(y_,ne0_)

  # plt.figure(41,figsize=(h_size, v_size))
  # plt.plot(y_,dphi_)
  
  #plt.figure(42,figsize=(h_size, v_size))
  #plt.plot(y_,vi1_)
  #plt.plot(y_,ve1_EXB_)
  #plt.plot(y_,ve1_D_)

  #plt.clear()
  plt.clf()
  legend_list=[]
  #plt.figure(43,figsize=(h_size, v_size))
  plt.plot(y_,Je_De_)
  legend_list = legend_list + ["Ie_De"] 
  plt.plot(y_,Ji0_)
  legend_list = legend_list + ["Ii0"] 
  plt.plot(y_, Ji1_)
  legend_list = legend_list + ["Ii1"]
  plt.plot(y_, JeExB_)
  legend_list = legend_list + ["IeExB"]
  plt.plot(y_, JeD_)
  legend_list = legend_list + ["IeD"]
  plt.plot(y_,Jtot_,'k',linewidth = linesize_large)
  legend_list = legend_list + ["Itot"]
  #plt.ylim((-3e-2, 3e-2))
  temp = integrate.cumtrapz(Jtot_, y_, initial = 0 )
  mean_current = temp[-1]
  plt.title( "{:.3f}".format(np.mean(Jtot_)*A_plasma*1000) + " mA" )
  #plt.title( "{:.3f}".format(mean_current*1000) + " mA" )
  plt.legend(legend_list, fancybox = False, edgecolor = 'k' ,loc  = 'right')


N_segments = len(m_decr)
segment_segments_ = 8
def animate_segments(ind):
  # ind = find_nearest_index(ion_velocity0,-500)
  m_decr_ = m_decr[ind]
  ky_ = m_decr_/R0
  f_  = frequency_decr[ind]
  T_ = 1/f_
  dt_ = m_decr_*T_/(N_segments-1)
  
  lambda_ = 2*np.pi/ky_
  segment_y = 2*np.pi*R0/segment_segments_
  
  y_ = np.linspace(0,segment_y,N_segments)
  y_2 = np.linspace(segment_y,2*segment_y,N_segments)
  dy_ = 2*np.pi*R0/(N_segments-1)
  
  De_ = electron_diffusion[ind]

  dphi_phase_ = dphi_phase_decr[ind]
  vi1_phase_ = vi1_phase_decr[ind]
  ve1_EXB_phase_ = ve1_EXB_phase_decr[ind] 
  ve1_D_phase_ = ve1_D_phase_decr[ind]
  
  t = np.linspace(0,T_ ,N_segments)
  Itot_v = np.zeros(N_segments)
  Itot_v2 = np.zeros(N_segments)
  
  # compute integral from 0 to pi/4
  i = 0
  t_ = 0

  
  while t_ <=  m_decr_*T_ + dt_/2:
    
    
    ne0_ = ne0[ind]*(1.0+np.sin(ky_*y_ - 2*np.pi*f_*t_))
    dphi_ = plasma_potential[ind] + dphi_decr[ind] * np.sin(ky_*y_ + dphi_phase_ - 2*np.pi*f_*t_)
    vi0_ = ion_velocity0[ind]
    vi1_ = vi1_decr[ind]*np.sin(ky_*y_ + vi1_phase_ - 2*np.pi*f_*t_)
    ve1_EXB_ = ve1_EXB_decr[ind]*np.sin(ky_*y_ + ve1_EXB_phase_ - 2*np.pi*f_*t_)
    ve1_D_ = ve1_D_decr[ind]*np.sin(ky_*y_ + ve1_D_phase_ - 2*np.pi*f_*t_)
    ve_De_ = De_/Lz  + kz*De_*np.sin(ky_*y_ + np.pi/2 - 2*np.pi*f_*t_)

    Je_De_ = -e*ne0_*ve_De_
    Ji0_ = e*ne0_*vi0_
    Ji1_ =  e*ne0_*vi1_
    JeExB_ = -e*ne0_*ve1_EXB_
    JeD_ = -e*ne0_*ve1_D_
    Jtot_1 = Ji0_ + Ji1_ + JeExB_ + JeD_ + Je_De_

    ne0_ = ne0[ind]*(1.0+np.sin(ky_*y_2 - 2*np.pi*f_*t_))
    dphi_ = plasma_potential[ind] + dphi_decr[ind] * np.sin(ky_*y_2 + dphi_phase_ - 2*np.pi*f_*t_)
    vi1_ = vi1_decr[ind]*np.sin(ky_*y_2 + vi1_phase_ - 2*np.pi*f_*t_)
    ve1_EXB_ = ve1_EXB_decr[ind]*np.sin(ky_*y_2 + ve1_EXB_phase_ - 2*np.pi*f_*t_)
    ve1_D_ = ve1_D_decr[ind]*np.sin(ky_*y_2 + ve1_D_phase_ - 2*np.pi*f_*t_)
    ve_De_ = De_/Lz  + kz*De_*np.sin(ky_*y_2 + np.pi/2 - 2*np.pi*f_*t_)

    Je_De_ = -e*ne0_*ve_De_
    Ji0_ = e*ne0_*vi0_
    Ji1_ =  e*ne0_*vi1_
    JeExB_ = -e*ne0_*ve1_EXB_
    JeD_ = -e*ne0_*ve1_D_
    Jtot_2 = Ji0_ + Ji1_ + JeExB_ + JeD_ + Je_De_

    #input("Enter any key")

    temp = integrate.cumtrapz(Jtot_1, y_, initial = 0 )
    Itot_v[i] = temp[-1]- temp[0]

    temp = integrate.cumtrapz(Jtot_2, y_2, initial = 0 )
    Itot_v2[i] = temp[-1]- temp[0]

    
    #print(temp[-1])
    i += 1
    t_ += dt_
  

  # plt.figure(40,figsize=(h_size, v_size))
  # plt.plot(y_,ne0_)

  # plt.figure(41,figsize=(h_size, v_size))
  # plt.plot(y_,dphi_)
  
  #plt.figure(44,figsize=(h_size, v_size))
  #plt.plot(y_,vi1_)
  #plt.plot(y_,ve1_EXB_)
  #plt.plot(y_,ve1_D_)

  plt.clf()
  legend_list=[]
  plt.plot(t, Itot_v,'b',linewidth = linesize_large)
  legend_list = legend_list + ["seg 1"]
  plt.plot(t, Itot_v2,'r',linewidth = linesize_large)
  legend_list = legend_list + ["seg 2"]
  plt.title( "{:.3f}".format(np.mean(Itot_v)*8) + " mA" )
  #plt.ylim((-4e-2, 4e-2))
  
  plt.legend(legend_list, fancybox = False, edgecolor = 'k' ,loc  = 'right')

if do_movies:
  fig = plt.figure(43,figsize=(h_size, v_size))
  ani = animation.FuncAnimation(fig, animate_currents, frames = range(0,len(m_decr)),interval = 80, repeat_delay=20)
  print(movies_dir + "movie_current.mp4")
  ani.save(movies_dir + "movie_current.mp4")
  fig = plt.figure(44,figsize=(h_size, v_size))
  ani = animation.FuncAnimation(fig, animate_segments, frames = range(0,len(m_decr)),interval = 80, repeat_delay=20)
  ani.save(movies_dir + "movie_segments.mp4")



ind = find_nearest_index(ion_velocity0,-900)
fig = plt.figure(20,figsize=(h_size, v_size))
animate_currents(ind)
plt.savefig(figures_dir+"currents_y.png", additional_artists = art, bbox_inches="tight")
plt.savefig(figures_dir+"currents_y.eps", additional_artists = art, bbox_inches="tight")

fig = plt.figure(21,figsize=(h_size, v_size))
animate_segments(ind)
plt.savefig(figures_dir+"segments_y.png", additional_artists = art, bbox_inches="tight")
plt.savefig(figures_dir+"segments_y.eps", additional_artists = art, bbox_inches="tight")


#                 W                       #


with open("linear_model_solutions/w_solution.txt") as s:
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
  
  plt.savefig(figures_dir+"wR_vs_ky.png")
  plt.savefig(figures_dir+"wR_vs_ky.eps")
  
  plt.figure(8,figsize=(h_size, v_size))
  

  plt.plot(ky*R0,wI*wI, linewidth = linesize)
  plt.plot(ky*R0,alpha, linewidth = linesize)
  plt.plot(ky*R0,beta, linewidth = linesize)
  plt.plot(ky*R0,gamma, linewidth = linesize)
  plt.legend(["$\omega{_I}^2$", "$α$", "$β$","$γ$"], loc = 'upper right' , fancybox = False, edgecolor = 'k' )
  plt.xlabel("m = k$_y$R$_0$")
  plt.ylabel("[rad/s]$^2$")
  plt.savefig(figures_dir+"wI_vs_ky.png")
  plt.savefig(figures_dir+"wI_vs_ky.eps")
  
  plt.figure(9,figsize=(h_size, v_size))
  plt.plot(ky*R0,[alpha[i]+beta[i]+gamma[i] for i in range(len(alpha))],'r', linewidth = linesize_large)
  plt.plot(ky*R0,[alpha[i]-beta[i]-gamma[i] for i in range(len(alpha))],'b', linewidth = linesize_large)
  plt.legend(["C$^2$", "$\omega{_I}^2$"], loc = 'upper right', fancybox = False, edgecolor = 'k' )
  plt.xlabel("m = k$_y$R$_0$")
  plt.ylabel("[rad/s]$^2$")
  plt.tight_layout()
  
  plt.savefig(figures_dir+"sqrt_vs_ky.png")
  plt.savefig(figures_dir+"sqrt_vs_ky.eps")
  
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
  
  plt.savefig(figures_dir+"sqrt_vs_ky.png")
  plt.savefig(figures_dir+"sqrt_vs_ky.eps")
  
plt.show()

