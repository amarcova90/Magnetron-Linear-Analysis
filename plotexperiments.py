from scipy import signal
import numpy as np
from numpy import unravel_index
import matplotlib.pyplot as plt
import csv
import glob
import os
import sys
import pandas as pd
from scipy.fftpack import fft
import pywt

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
do_show = False
do_save = True
do_spectrum = False

date = "2021-06-18"
file_pramble = "Y" #"WFM"
run_all = False
#runs_subset = [1]
runs_subset = range(1,33+1)
excelfile = "experiment/Testing_2020.xlsx"

df = pd.read_excel(excelfile, sheet_name = date)

# values to save (one for run)
run_list = []
I1 = []
I2 = []
I3 = []
I_tot = []
dV_ballast = []
mode = []
main_frequency = []


# start routine
for row,seg1 in enumerate(df["file name current 1"]):
  if (type(seg1) is int) and (df["Run"][row] in runs_subset):
    seg1 = str(seg1)
    seg2 = str(df["file name current 2"][row])
    seg3 = str(df["file name current 3"][row])
    seg4 = str(df["file name current 4"][row])

    gas_name = df["gas"][row]
    data_location = "experiment/"+gas_name+"/"+date
    run_num = str(df["Run"][row])
    if len(seg1) == 1:
      seg1 = "0" + seg1
    if len(seg2) == 1:
      seg2 = "0" + seg2
    if len(seg3) == 1:
      seg3 = "0" + seg3
    if len(seg4) == 1:
      seg4 = "0" + seg4      
    s1_file = data_location + "/" + file_pramble + seg1 + ".csv"
    s2_file = data_location + "/" + file_pramble + seg2 + ".csv"
    s3_file = data_location + "/" + file_pramble + seg3 + ".csv"
    s4_file = data_location + "/" + file_pramble + seg4 + ".csv"

    # print(s1_file)
    # print(s2_file)
    # print(s3_file)
    # print(s4_file)
    # print("\n")
    # input()

    t = []
    s1 = []
    s2 = []
    s3 = []
    s4 = []

    with open(s1_file) as csvfile:
      csvreader = csv.reader(csvfile)
      next(csvreader)
      for row in csvreader:
        #print(row)
        t.append(float(row[0]))
        s1.append(float(row[1]))

    with open(s2_file) as csvfile:
      csvreader = csv.reader(csvfile)
      next(csvreader)
      for row in csvreader:
        s2.append(float(row[1]))

    with open(s3_file) as csvfile:
      csvreader = csv.reader(csvfile)
      next(csvreader)
      for row in csvreader:
        s3.append(float(row[1]))

    with open(s4_file) as csvfile:
      csvreader = csv.reader(csvfile)
      next(csvreader)
      for row in csvreader:
        s4.append(float(row[1]))

    t = np.array(t)
    s1 = np.array(s1)
    s2 = np.array(s2)
    s3 = np.array(s3)
    s4 = np.array(s4)
    
    Itot = s1 + s2 + s3
    
    Ln = len(s1)

    
    
    Nma = 25 # only odd
    nma = int(Nma/2)
    s1_ma = np.zeros(Ln-2*nma)
    s2_ma = np.zeros(Ln-2*nma)
    s3_ma = np.zeros(Ln-2*nma)
    t_ma = t[nma:-(nma)] 
    
    k = 0
    for i in range(nma,Ln-nma):
      s1_ma[k] = s1[i]
      s2_ma[k] = s2[i]
      s3_ma[k] = s3[i]
      for j in range(1,nma+1):
        s1_ma[k] = s1_ma[k] + s1[i-j] + s1[i+j]
        s2_ma[k] = s2_ma[k] + s2[i-j] + s2[i+j]
        s3_ma[k] = s3_ma[k] + s3[i-j] + s3[i+j]
      s1_ma[k] = s1_ma[k]/Nma
      s2_ma[k] = s2_ma[k]/Nma
      s3_ma[k] = s3_ma[k]/Nma
      k = k + 1
      
    Itot_ma = s1_ma + s2_ma + s3_ma

    print("ballast voltage drop = {}".format(np.mean(s4)*500)) 

    # remove mean value
    # s1 = s1 - np.mean(s1)
    # s2 = s2 - np.mean(s2)

    # normalize
    #s1 = s1/np.max(s1)
    #s2 = s2/np.max(s2)
    # s1 = s1/np.max(s1)

    s1_detrended = signal.detrend(s1)
    s2_detrended = signal.detrend(s2)
    s3_detrended = signal.detrend(s3)

    R0 = 9.5e-3
    #dx = 10.9*np.pi/180.0
    dx = 10*np.pi/180.0
    # dx = 20*np.pi/180.0
    #dx = 40.0*np.pi/180.0

    dt = t[1] - t[0]
    print(dt)
    #input()
    # wk2p01(s1,s2,dt,dx)


    ff1 = fft(s1_detrended)
    ff2 = fft(s2_detrended)
    ff1 = ff1[0:len(ff1)//2]
    ff2 = ff2[0:len(ff2)//2]


    fs = 1.0/dt
    f_v = fs*np.arange(0.0,len(ff1),1.0)/len(ff1)/2.0

    max_f = np.amax(np.abs(ff1))
    indx = np.where(np.abs(ff1) == max_f)


    w12 = ff1*np.conj(ff2)
    #print(w12)
    kx12 = np.angle(w12)/dx

    # print("main mode = {}".format(kx12[indx]))
    # print("main frequency = {}".format(f_v[indx]))

    # plt.figure(figsize=(h_size, v_size))
    # plt.plot(f_v,kx12)
      
    ## spectrum analysis
    if do_spectrum:

      nf = 1000
      nk = 1000
      n = len(s1_detrended)

      kNyq = 2*np.pi/dx/2.0			# Nyquist wavenumber
      fNyq = 1.0/dt/2.0		# Nyquist frequency
      fupper = fNyq/2.001	# maximum resolvable frequency
      flower = 4.0/dt/n		# minimum resolvable frequency
      #flower = 1e6

      f = np.logspace(np.log10(flower),np.log10(fupper),nf)

      #f = np.linspace(3e6,6e6,nf)

      #f = np.linspace(flower,fupper,nf)
      #print(f)
      scale = 1/(f*dt)
      fz = pywt.scale2frequency('morl', scale)/dt
      k = np.linspace(-kNyq,kNyq,nk)
      dk = 2*kNyq / nk
      S = np.zeros((nf,nk))
      H = np.zeros((nf,nk))
      kur = np.zeros((nf,nk))
      S11 = np.zeros((nf,nk),dtype = complex)
      S12 = np.zeros((nf,nk),dtype = complex)
      S22 = np.zeros((nf,nk),dtype = complex)
      fsp = np.zeros(nf)

      wx1 = signal.cwt(s1_detrended, signal.morlet2, scale) 

      wx2 = signal.cwt(s2_detrended, signal.morlet2, scale)

      print("len = {}".format(len(wx1)))

      for i in range(nf):
        print("spectral analysis {:0.1f} %".format(i/nf*100))
        ntrans = int(round(np.min([1.5*scale[i], n/8])))
        #print(ntrans)
        
        wwx1 = wx1[i][ntrans:-ntrans] # remove transients
        wwx2 = wx2[i][ntrans:-ntrans] # remove transients

        w11 = wwx1*np.conj(wwx1)
        w22 = wwx2*np.conj(wwx2)
        w12 = wwx1*np.conj(wwx2)
        kx = np.angle(w12)/dx
        #print(kx)
        neff = len(kx)
        #print(neff[i])
        ampl = np.abs((w11+w22) / 2.0 /float(neff))
        ampl2 = ampl*ampl
        pos = np.abs((kx+kNyq))/(dk+sys.float_info.epsilon)
        #print(kx)
        fsp[i] = np.sum(ampl) / neff
        for j in range(neff):
          ind = int(pos[j])
          H[i,ind] = H[i,ind] + 1
          S[i,ind] = S[i,ind] + ampl[j]
          kur[i,ind] = kur[i,ind] + ampl2[j]
          S11[i,ind] = S11[i,ind] + w11[j]
          S22[i,ind] = S22[i,ind] + w22[j]
          S12[i,ind] = S12[i,ind] + w12[j]
          
          

      i_f, i_t = np.unravel_index(np.argmax(S), S.shape)
      m_main = k[i_t]
      max_f = fz[i_f]

      plt.figure(figsize=(h_size, v_size))
      plt.pcolormesh(k,fz*1e-3,np.log10(S), cmap = 'jet')
      #plt.plot(m_main,max_f*1e-3,'ko')
      #plt.pcolormesh(k,fz*1e-3,S, cmap = 'jet')
      # plt.pcolormesh(k + 2*kNyq,fz*1e-3,np.log10(S), cmap = 'jet')
      # plt.pcolormesh(k - 2*kNyq,fz*1e-3,np.log10(S), cmap = 'jet')
      plt.xlabel("mode")
      plt.ylabel("frequency [kHz]")
      plt.ylim((0,10*1e3))
      plt.colorbar()
      plt.tight_layout()
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_spectrum.png") 
        plt.close()
      

      plt.figure(figsize=(h_size, v_size))
      plt.plot(f_v*1e-3,np.abs(ff1))
      plt.plot(f_v*1e-3,np.abs(ff2))
      plt.xlabel("frequency [kHz]")
      plt.ylabel("signal amp [mA]")
      plt.legend(["s1","s2"])
      plt.tight_layout()
      plt.xlim((0,10*1e3))
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_fft.png") 
        plt.close()
      # filename = os.path.join(path,"FFT.png") 
      # plt.savefig(filename)  

      plt.figure(figsize=(h_size, v_size))
      plt.subplot(2,1,1)
      plt.semilogy(f_v/1e6,np.abs(ff1),'b')
      plt.tight_layout() 
      plt.axis((0,10,1e-4,1))
      plt.subplot(2,1,2)
      plt.semilogy(f_v/1e6,np.abs(ff2),'r')
      plt.tight_layout() 
      plt.axis((0,10,1e-4,1))
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_fft_semilogy.png") 
        plt.close()
 
      nf = len(ff1)

      f = np.linspace(flower,fupper,nf) # array of frequencies
      scale = 1.0/(f*dt) # array of wavelet scales

      Z1 = signal.cwt(s1_detrended, signal.morlet2, scale) 

      Z2 = signal.cwt(s2_detrended, signal.morlet2, scale) 

      fz = pywt.scale2frequency('morl', scale)/dt


      plt.figure(figsize=(h_size, v_size))
      plt.contourf(t*1e6,fz*1e-3,np.abs(Z1),levels = 100, cmap = 'jet')
      plt.colorbar()
      plt.ylim((0,10*1e3))
      plt.xlabel("time [$\mu$s]")
      plt.ylabel("frequency [kHz]")
      plt.title("signal 1")
      plt.tight_layout()
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_wavelet_s1.png") 
        plt.close()
      # filename = os.path.join(path,"signal1Morlet_contour.png") 
      # plt.savefig(filename) 

      i_f, i_t = np.unravel_index(np.argmax(Z1), Z1.shape)



      plt.figure(figsize=(h_size, v_size))
      plt.contourf(t*1e6,fz*1e-3,np.abs(Z2),levels = 100, cmap = 'jet')
      plt.colorbar()
      plt.ylim((0,10*1e3))
      plt.xlabel("time [$\mu$s]")
      plt.ylabel("frequency [kHz]")
      plt.title("signal 2")
      plt.tight_layout() 
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_wavelet_s2.png") 
        plt.close()
      # filename = os.path.join(path,"signal2Morlet_contour.png") 
      # plt.savefig(filename) 
 
      
    fig,ax = plt.subplots(figsize=(h_size, v_size))
    lns1 = ax.plot(t_ma*1e6,s1_ma/50*1000, label = '$s_1$')
    lns2 = ax.plot(t_ma*1e6,s2_ma/50*1000, label = '$s_2$')
    ax.set_xlabel("time [$\mu$s]")
    ax.set_ylabel("signal current 1,2 [mA]")
    ax.set_ylim((np.mean(s2/50*1000)*5, np.max([np.max(s1), np.max(s2)])))
    ax2=ax.twinx()
    lns3  = ax2.plot(t_ma*1e6,s3_ma/50*1000, label = '$s_3$')
    lns4  = ax2.plot(t_ma*1e6,Itot_ma/50*1000, "--k", label = '$I_{tot}$')
    lns = lns1 + lns2 + lns3 + lns4
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs)#, loc=7)
    #ax2.set_ylim((np.mean(s3/50*1000) + (np.mean(s1/50*1000) + (np.mean(s1/50*1000)))*4,0.1))
    ax2.set_ylim((np.mean(Itot_ma/50*1000)*1.05,0.1))
    ax2.set_ylabel("signal current 3 (main), total [mA]")
    plt.tight_layout()
    plt.title("signals")
    if do_save:
      plt.savefig(data_location+"/"+run_num+"_signal_filtered.png") 
      plt.close()
    

    plt.figure(figsize=(h_size, v_size))
    plt.plot(t*1e6,s1/50*1000)
    plt.plot(t*1e6,s2/50*1000)
    plt.xlabel("time [$\mu$s]")
    plt.ylabel("signal current [mA]")
    plt.legend(["s1","s2"])
    #plt.xlim((0,10))
    plt.title("raw signals")
    if do_save:
      plt.savefig(data_location+"/"+run_num+"_signal_raw_12.png") 
      plt.close()



    fig,ax = plt.subplots(figsize=(h_size, v_size))
    lns1 = ax.plot(t*1e6,s1/50*1000, label = 's1')
    lns2 = ax.plot(t*1e6,s2/50*1000, label = 's2')
    ax.set_xlabel("time [$\mu$s]")
    ax.set_ylabel("signal current 1,2 [mA]")
    ax.set_ylim((-0.65,0.05))
    ax2=ax.twinx()
    lns3  = ax2.plot(t*1e6,s3/50*1000, "k", label = 's3')
    lns = lns1 + lns2 + lns3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=7)
    ax2.set_ylim((np.mean(s3/50*1000)-0.5,0.1))
    ax2.set_ylabel("signal current 3 (main) [mA]")
    plt.tight_layout()
    plt.title("raw signals")
    if do_save:
      plt.savefig(data_location+"/"+run_num+"_signal_raw_all.png")
      plt.close()
      




    
    
    print("Currents 1, 2, 3, sum [mA]:")
    print(np.mean(s1)/50*1000)
    print(np.mean(s2)/50*1000)
    print(np.mean(s3)/50*1000)
    print(np.mean(s1)/50*1000+np.mean(s2)/50*1000+np.mean(s3)/50*1000)
    print("ballast voltage drop = {}".format(np.mean(s4)*500))  

    if do_spectrum:
      print("Max frequency, mode:")
      print("main mode = {}".format(m_main))
      print("main frequency = {}".format(max_f))   
    
    # add to lists
    run_list.append(run_num)
    I1.append(-np.mean(s1)/50*1000)
    I2.append(-np.mean(s2)/50*1000)
    I3.append(-np.mean(s3)/50*1000)
    I_tot.append(-np.mean(s1)/50*1000-np.mean(s2)/50*1000-np.mean(s3)/50*1000)
    dV_ballast.append(np.mean(s4)*500)
    if do_spectrum:
      mode.append(m_main)
      main_frequency.append(max_f)




print("run #:")
print(run_list)
print("I1 mean:")
print(I1)
print("I2 mean:")
print(I2)
print("I3 mean:")
print(I3)
print("Itot mean:")
print(I_tot)
print("dV ballast:")
print(dV_ballast)
print("main mode:")
print(mode)
print("main frequency:")
print(main_frequency)

# f = open(data_location + "/solution.txt", "w")

# f.writelines(run_list)
# f.writelines(map(str, I1))
# f.writelines(map(str, I2))
# f.writelines(map(str, I3))
# f.writelines(map(str, I_tot))
# f.writelines(map(str, dV_ballast))
# f.writelines(map(str, mode))
# f.writelines(map(str, main_frequency))

# f.close()

if do_show:
  plt.show()
