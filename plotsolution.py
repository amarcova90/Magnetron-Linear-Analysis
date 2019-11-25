import matplotlib.pyplot as plt
import numpy as np

with open("solution.txt") as s:
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
				i = 0
			else:
				try:
					temp_value = float(temp[0])
					plasma_potential[i] = temp_value
					m_decr[i] = float(temp[1])
					frequency_decr[i] = float(temp[2])
					Ix_decr[i] = float(temp[3])
					m_incr[i] = float(temp[4])
					frequency_incr[i] = float(temp[5])
					Ix_incr[i] = float(temp[6])
					
					i = i + 1
				except(ValueError):
					pass
		except(IndexError):
			pass

unique_m = set(list(set(m_decr)) + list(set(m_incr)))
legend_list=[]
colors = {1: (0.6350, 0.0780, 0.1840),2: (0, 0.4470, 0.7410),3: (0.8500, 0.3250, 0.0980), 3: (0.9290, 0.6940, 0.1250), 4: (0.4940, 0.1840, 0.5560), 5: (0.4660, 0.6740, 0.1880), 6 : (0.3010, 0.7450, 0.9330), 7 : (1 , 0, 0)}

plt.figure(1)
for mode in unique_m:
	if mode in m_decr:
		plt.plot(abs(plasma_potential[m_decr == mode]),frequency_decr[m_decr == mode]/1000,color=colors[mode])
		legend_list = legend_list + ["m$_{decr}$ = " + str(mode)] 
	if mode in m_incr:
		plt.plot(abs(plasma_potential[m_incr == mode]),frequency_incr[m_incr == mode]/1000,'--',color = colors[mode])	
		legend_list = legend_list + ["m$_{incr}$ = " + str(mode)]
		
plt.legend(legend_list)
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Frequency [kHz]")
#plt.show(block=True)
plt.show()
plt.savefig("f_vs_EL_both.png")


plt.figure(2)
for mode in unique_m:
	if mode in m_decr:
		plt.plot(abs(plasma_potential[m_decr == mode]),Ix_decr[m_decr == mode]*1000,color=colors[mode])
		legend_list = legend_list + ["m_decr = " + str(mode)] 
	if mode in m_incr:
		plt.plot(abs(plasma_potential[m_incr == mode]),Ix_incr[m_incr == mode]*1000,'--',color = colors[mode])	
		legend_list = legend_list + ["m_incr = " + str(mode)]
		
plt.legend(legend_list)
plt.xlabel("E$_0$L$_n$ [V]")
plt.ylabel("Current [mA]")
#plt.show(block=True)
plt.show()
plt.savefig("I_vs_EL_both.png")

	

    
