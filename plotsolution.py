import matplotlib.pyplot as plt
import numpy as np
import os
import sys

if len(sys.argv)==2:
	solution_file=sys.argv[1]
	isexperiment=0
elif len(sys.argv)==3:
	solution_file=sys.argv[1]
	data_file=sys.argv[2]
	isexperiment=1
else:
	print("Wrong number of inputs. Usage: python3 postprocess.py <solution_file> <data_file> (optional)")
	sys.exit()



with open(solution_file) as s:
	lines = s.read().splitlines()
	inputs = ("ion_mass_u", "Te", "B0", "R0", "LB", "Ln", "kz", "kx", "N_voltages")
	for line in lines:
		# first I read the input parameters
		# I need to know N_voltages to initialize
		# the numpy arrays.
		temp=line.split()
		try:
			#print(temp[1])
			if temp[0] == 'ion_mass_u':
				ion_mass_u = float(temp[1])
			elif temp[0] == "Te":
				Te = float(temp[1])
			elif temp[0] == "B0":
				B0 = float(temp[1])
			elif temp[0] == "R0":
				R0 = float(temp[1])
			elif temp[0] == "LB":
				LB = float(temp[1])
			elif temp[0] == "Ln":
				Ln = float(temp[1])
			elif temp[0] == "kz":
				kz = float(temp[1])
			elif temp[0] == "kx":
				kx = float(temp[1])
			elif temp[0] == "N_voltages":
				N_voltages = int(temp[1])
				plasma_potential = np.zeros(N_voltages, dtype = float)
				m_decr = np.zeros(N_voltages, dtype = int)	
				frequency_decr = np.zeros(N_voltages, dtype = float)
				Ix_decr = np.zeros(N_voltages, dtype = float)
				m_incr = np.zeros(N_voltages, dtype = int)					
				frequency_incr = np.zeros(N_voltages, dtype = float)
				Ix_incr = np.zeros(N_voltages, dtype = float)
				ne_over_n0_decr = np.zeros(N_voltages, dtype = float)
				ne_over_n0_incr = np.zeros(N_voltages, dtype = float)
				vi1_decr = np.zeros(N_voltages, dtype = float)
				vi1_incr = np.zeros(N_voltages, dtype = float)
				ion_velocity0 = np.zeros(N_voltages, dtype = float)
				i = 0
			else:
				try:
					temp_value = float(temp[0])
					plasma_potential[i] = temp_value
					m_decr[i] = float(temp[1])
					frequency_decr[i] = float(temp[2])
					vi1_decr[i] = float(temp[3])
					Ix_decr[i] = float(temp[4])
					ne_over_n0_decr[i] = float(temp[5])
					m_incr[i] = float(temp[6])
					frequency_incr[i] = float(temp[7])
					vi1_incr[i] = float(temp[8])
					Ix_incr[i] = float(temp[9])
					ne_over_n0_incr[i] = float(temp[10])
					ion_velocity0[i] = float(temp[11])
					
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
				gap_voltage_exp[i-1]=float(temp[0])
				gap_current_exp[i-1]=float(temp[1])
				mode_number_exp[i-1]=float(temp[2])
				frequency_exp[i-1]=float(temp[3])
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
		voltage_shift = float(lines[-1].split()[0])
		
			
unique_m = set(list(set(m_decr)) + list(set(m_incr)))
legend_list=[]
colors = {np.nan: 'k', 1: (0.6350, 0.0780, 0.1840),2: (0, 0.4470, 0.7410),3: (0.8500, 0.3250, 0.0980), 3: (0.9290, 0.6940, 0.1250), 4: (0.4940, 0.1840, 0.5560), 5: (0.4660, 0.6740, 0.1880), 6 : (0.3010, 0.7450, 0.9330), 7 : (1 , 0, 0)}

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, solution_file[:-4] + '_figures/')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
	
plt.figure(1)
for mode in unique_m:
	if mode in m_decr:
		plt.plot(abs(plasma_potential[m_decr == mode]),frequency_decr[m_decr == mode]/1000,color=colors[mode])
		legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
	if mode in m_incr:
		plt.plot(abs(plasma_potential[m_incr == mode]),frequency_incr[m_incr == mode]/1000,'--',color = colors[mode])	
		legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]

if isexperiment:
	for mode in set(mode_number_exp_decr[~np.isnan(mode_number_exp_decr)]):
		plt.plot(abs(gap_voltage_exp_decr[mode_number_exp_decr == mode])-voltage_shift,frequency_exp_decr[mode_number_exp_decr == mode]/1000,'v',color = colors[mode])	
		legend_list = legend_list + ["m$_{decr}$ = " + str(int(mode))]
	for mode in set(mode_number_exp_incr[~np.isnan(mode_number_exp_incr)]):
		plt.plot(abs(gap_voltage_exp_incr[mode_number_exp_incr == mode])-voltage_shift,frequency_exp_incr[mode_number_exp_incr == mode]/1000,'^',color = colors[mode])	
		legend_list = legend_list + ["m$_{incr}$ = " + str(int(mode))]

		
		
plt.legend(legend_list,loc="upper right" )
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Frequency [kHz]")
#plt.show(block=True)
plt.show()
plt.savefig(results_dir+"f_vs_EL.png")


plt.figure(2)
for mode in unique_m:
	if mode in m_decr:
		plt.plot(abs(plasma_potential[m_decr == mode]),Ix_decr[m_decr == mode]*1000,color=colors[mode])
		legend_list = legend_list + ["m_{decr} = " + str(mode)] 
	if mode in m_incr:
		plt.plot(abs(plasma_potential[m_incr == mode]),Ix_incr[m_incr == mode]*1000,'--',color = colors[mode])	
		legend_list = legend_list + ["m_{incr} = " + str(mode)]

if isexperiment:
	for mode in set(mode_number_exp_decr[~np.isnan(mode_number_exp_decr)]):
		plt.plot(abs(gap_voltage_exp_decr[mode_number_exp_decr == mode])-voltage_shift,gap_current_exp_decr[mode_number_exp_decr == mode],'v',color = colors[mode])	
		legend_list = legend_list + ["m$_{decr}$ = " + str(int(mode))]
	for mode in set(mode_number_exp_incr[~np.isnan(mode_number_exp_incr)]):
		plt.plot(abs(gap_voltage_exp_incr[mode_number_exp_incr == mode])-voltage_shift,gap_current_exp_incr[mode_number_exp_incr == mode],'^',color = colors[mode])	
		legend_list = legend_list + ["m$_{incr}$ = " + str(int(mode))]

plt.legend(legend_list)
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Current [mA]")
#plt.show(block=True)
plt.show()
plt.savefig(results_dir+"I_vs_EL.png")


plt.figure(3)
for mode in unique_m:
	if mode in m_decr:
		plt.plot(abs(plasma_potential[m_decr == mode]),ne_over_n0_decr[m_decr == mode],color=colors[mode])
		legend_list = legend_list + ["m_{decr} = " + str(mode)] 
	if mode in m_incr:
		plt.plot(abs(plasma_potential[m_incr == mode]),ne_over_n0_incr[m_incr == mode],'--',color = colors[mode])	
		legend_list = legend_list + ["m_{incr} = " + str(mode)]

plt.legend(legend_list)
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("ne/n0 []")
#plt.show(block=True)
plt.show()
plt.savefig(results_dir+"ne_over_n0_vs_EL.png")

plt.figure(4)
for mode in unique_m:
	if mode in m_decr:
		plt.plot(abs(plasma_potential[m_decr == mode]),-vi1_decr[m_decr == mode],color=colors[mode])
		legend_list = legend_list + ["m_{decr} = " + str(mode)] 
	if mode in m_incr:
		plt.plot(abs(plasma_potential[m_incr == mode]),-vi1_incr[m_incr == mode],'--',color = colors[mode])	
		legend_list = legend_list + ["m_{incr} = " + str(mode)]

plt.legend(legend_list)
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("v$_{i1}$ [m/s]")
#plt.show(block=True)
plt.show()
plt.savefig(results_dir+"vi1_vs_EL.png")
    
