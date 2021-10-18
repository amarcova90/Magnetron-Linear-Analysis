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
from matplotlib.ticker import FormatStrFormatter

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
do_spectrum = True

nf = 1000
nk = 1000

f_ceil = 10*1e3
date = "2021-02-21"
file_pramble = "WF_P"
run_all = False
#runs_subset = [58]
runs_subset = range(1,150+1)
excelfile = "experiment/Testing_2020_new3.xlsx"

df = pd.read_excel(excelfile, sheet_name = date)

# values to save (one for run)
run_list = []
I1 = []
I2 = []
I3 = []
I_tot = []
I_tot_mean = []
V_cathode = []
V_cathode_mean = []
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
    t4 = []
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
        if float(float(row[1])<0.0):
          s4.append(-float(row[1]))
        else:
          s4.append(float(row[1]))
        t4.append(float(row[0]))

    t = np.array(t)
    s1 = np.array(s1)
    s2 = np.array(s2)
    s3 = np.array(s3)
    s4 = np.array(s4)
    
    Itot = (s1 + s2 + s3)/50.0
    
    Ln = len(s1)

   # mean average filter 
    
    Nma = 25 # only odd numbers
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


    # remove mean value
    # s1 = s1 - np.mean(s1)
    # s2 = s2 - np.mean(s2)

    # normalize
    #s1 = s1/np.max(s1)
    #s2 = s2/np.max(s2)
    # s1 = s1/np.max(s1)

    s1_detrended = signal.detrend(s1)
    s2_detrended = signal.detrend(s2)

    s1_mod = s1_detrended/np.max(np.abs(s1_detrended))
    s2_mod = s2_detrended/np.max(np.abs(s2_detrended))

    #plt.figure()
    #plt.plot(t, s1_mod)
    #plt.plot(t, s2_mod)
    #plt.show()
    #input()

    R0 = 9.5e-3
    #dx = 10.9*np.pi/180.0
    dx = 10*np.pi/180.0
    # dx = 20*np.pi/180.0
    #dx = 40.0*np.pi/180.0

    dt = t[1] - t[0]
    #print(dt)
    #input()
    # wk2p01(s1,s2,dt,dx)


    ff1 = fft(s1_mod)
    ff2 = fft(s2_mod)
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

      n = len(s1_mod)

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

      wx1 = signal.cwt(s1_mod, signal.morlet2, scale) 

      wx2 = signal.cwt(s2_mod, signal.morlet2, scale)

      for i in range(nf):
        print("Run {}: spectral analysis {:0.1f} %".format(run_num, i/nf*100))
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
      f_main = fz[i_f]    
      
      # find time of maximum 
      
#      plt.figure()
#      plt.plot(t, np.abs(wx1[i_f,:]))
      index_t1 = np.argmax(wx1[i_f,:])
      max_wx1 = np.abs(wx1[i_f,:])[index_t1]
#      plt.plot(t[index_t1], max_wx1, "ok")

#      plt.figure()
#      plt.plot(t, np.abs(wx2[i_f,:]))
      index_t2 = np.argmax(wx2[i_f,:])
      max_wx2 = np.abs(wx2[i_f,:])[index_t2]
#      plt.plot(t[index_t2], max_wx2, "ok")
     
      if np.abs(wx1[i_f,index_t1]) > np.abs(wx2[i_f,index_t2]):
        Vc_main = s4[index_t1]
        Itot_main = Itot[index_t1]
      else:
        Vc_main = s4[index_t2]
        Itot_main = Itot[index_t2]

      # plot maximum of spectrum at each frequency with same color scale as spectrum
      f_ = []
      s_ = []
      k_ = []

      for i in range(len(k)):
        k_.append(k[i])
        i_ = np.argmax(S[:,i])
        f_.append(fz[i_])
        s_.append(S[i_,i])
        #print(min(S[:,i]))

      f_.append(-2000000.0)
      k_.append(0.0)
      s_.append(100.0)

      f_ = np.array(f_)
      k_ = np.array(k_)

      s_ = np.log10(s_)
      s_[-1] = np.min(np.log10(S[~np.isinf(np.log10(S))]))
      


      fig = plt.figure(figsize=(h_size, v_size))
      plt.scatter(k_,f_*1e-3, c = s_, cmap = "jet")
      plt.xlabel("mode")
      plt.ylabel("frequency [kHz]")
      plt.ylim((0,10*1e3)) 
      plt.tight_layout()
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_maxf_m.png")
        fig.clear()
        plt.close()


      maxS = np.log10(np.max(S))
      fig = plt.figure(figsize=(h_size, v_size))
      plt.pcolormesh(k,fz*1e-3,np.log10(S)/maxS, cmap = 'jet', shading = 'auto')
      #plt.plot(m_main,f_main*1e-3,'bx', markersize  = 10)
      #plt.pcolormesh(k,fz*1e-3,S, cmap = 'jet')
      # plt.pcolormesh(k + 2*kNyq,fz*1e-3,np.log10(S), cmap = 'jet')
      # plt.pcolormesh(k - 2*kNyq,fz*1e-3,np.log10(S), cmap = 'jet')
      plt.xlabel("mode")
      plt.ylabel("frequency [kHz]")
      plt.ylim((0,f_ceil))
      plt.colorbar()
      plt.tight_layout()
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_spectrum.png") 
        fig.clear()
        plt.close()
      

      fig = plt.figure(figsize=(h_size, v_size))
      plt.plot(f_v*1e-3,np.abs(ff1))
      plt.plot(f_v*1e-3,np.abs(ff2))
      plt.xlabel("frequency [kHz]")
      plt.ylabel("signal amp [mA]")
      plt.legend(["s1","s2"])
      plt.tight_layout()
      plt.xlim((0,f_ceil))
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_fft.png") 
        fig.clear()
        plt.close()
      # filename = os.path.join(path,"FFT.png") 
      # plt.savefig(filename)  

      fig = plt.figure(figsize=(h_size, v_size))
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
        fig.clear()
        plt.close()

      maxwx = np.max(np.concatenate([np.abs(wx1),np.abs(wx2)]))
      #wx1[0,-1] = maxwx
      #wx2[0,-1] = maxwx

 #     plt.contourf(t*1e6,fz*1e-3,np.log10(np.abs(wx1)),levels = 100, cmap = 'jet')
      fig, ax = plt.subplots(figsize=(h_size, h_size))
      m = ax.pcolormesh(t*1e6,fz*1e-3,np.abs(wx1)/np.max(np.abs(wx1)), cmap = 'jet', shading ='auto')
      cbar = fig.colorbar(m)
      cbar.ax.tick_params(labelsize = 30)
      ax.set_ylim((0,f_ceil))
      ax.set_xlabel("time [$\mu$s]")
      ax.set_ylabel("frequency [kHz]")
      ax.set_title("segment 1")
      for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(30)
      plt.tight_layout()
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_wavelet_s1.png") 
        fig.clear()
        plt.close()
      # filename = os.path.join(path,"signal1Morlet_contour.png") 
      # plt.savefig(filename) 

      #i_f, i_t = np.unravel_index(np.argmax(Z1), Z1.shape)



      fig, ax = plt.subplots(figsize=(h_size, h_size))
      m = ax.pcolormesh(t*1e6,fz*1e-3,np.abs(wx2)/np.max(np.abs(wx2)), cmap = 'jet', shading = 'auto')
      cbar = fig.colorbar(m)
      cbar.ax.tick_params(labelsize = 30)
      ax.set_ylim((0,f_ceil))
      ax.set_xlabel("time [$\mu$s]")
      ax.set_ylabel("frequency [kHz]")
      ax.set_title("segment 2")
      for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
          item.set_fontsize(30)
      plt.tight_layout() 
      if do_save:
        plt.savefig(data_location+"/"+run_num+"_wavelet_s2.png") 
        fig.clear()
        plt.close()
      # filename = os.path.join(path,"signal2Morlet_contour.png") 
      # plt.savefig(filename) 



    ax1max = 0.0
    ax1min = np.min(np.concatenate([s1_ma,s2_ma]))*1.15/50*1000.0 
    ax2max = np.max(s3_ma)*0.67/50*1000.0
    ax2min = np.min(Itot_ma)/50*1000.0*1.08
      
    fig, (ax, ax2) = plt.subplots(2, 1, sharex=True, figsize=(h_size, v_size))
    lns1 = ax.plot(t_ma*1e6,s1_ma/50*1000, label = '$s_1$')
    lns2 = ax.plot(t_ma*1e6,s2_ma/50*1000, label = '$s_2$')
    
    ax.set_ylabel("segment current [mA]")
    #ax.set_ylim((np.mean(s2/50*1000)*6, np.max([np.max(s1), np.max(s2)])))
    #ax.set_ylim((np.mean(s2/50*1000)*6, 0.0))
    ax.set_ylim((ax1min, ax1max))
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    #ax2=ax.twinx()
    lns3  = ax2.plot(t_ma*1e6,s3_ma/50*1000, color = '#2ca02c', label = '$s_3$')
    lns4  = ax2.plot(t_ma*1e6,Itot_ma/50*1000, "-", color = '#d62728', label = '$total$')
    lns = lns1 + lns2 + lns3 + lns4
    labs = [l.get_label() for l in lns]
    
    #ax2.set_ylim((np.mean(s3/50*1000) + (np.mean(s1/50*1000) + (np.mean(s1/50*1000)))*4,0.1))
    #ax2.set_ylim((np.mean(Itot_ma/50*1000)*1.05,0.0))
    ax2.set_ylim((ax2min,ax2max))
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    ax2.set_xlabel("time [$\mu$s]")
    ax2.legend(lns, labs, ncol=2, loc=1)
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    #ax2.set_ylabel("signal current 3 (main), total [mA]")
    
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    
    plt.tight_layout()
    #plt.title("signals")
    if do_save:
      plt.savefig(data_location+"/"+run_num+"_signal_filtered.png") 
      fig.clear()
      plt.close()
    

    fig = plt.figure(figsize=(h_size, v_size))
    plt.plot(t*1e6,s1/50*1000)
    plt.plot(t*1e6,s2/50*1000)
    plt.xlabel("time [$\mu$s]")
    plt.ylabel("signal current [mA]")
    plt.legend(["s1","s2"])
    #plt.xlim((0,10))
    plt.title("raw signals")
    if do_save:
      plt.savefig(data_location+"/"+run_num+"_signal_raw_12.png") 
      fig.clear()
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
      fig.clear()
      plt.close()


    fig = plt.figure(figsize=(h_size, v_size))
    plt.plot(t*1e6,s4)
    plt.xlabel("time [$\mu$s]")
    plt.ylabel("voltage [mV]")
    plt.title("cathode voltage")
    if do_save:
      plt.savefig(data_location+"/"+run_num+"_cathode_voltage.png") 
      fig.clear()
      plt.close()

    fig = plt.figure(figsize=(h_size, v_size))
    plt.plot(t*1e6,Itot*1e3)
    plt.xlabel("time [$\mu$s]")
    plt.ylabel("current [mA]")
    plt.title("total current")
    if do_save:
      plt.savefig(data_location+"/"+run_num+"_total_current.png") 
      fig.clear()
      plt.close()

    
    print("Currents 1, 2, 3, sum [mA]:")
    print(np.mean(s1)/50*1000)
    print(np.mean(s2)/50*1000)
    print(np.mean(s3)/50*1000)
    print(np.mean(s1)/50*1000+np.mean(s2)/50*1000+np.mean(s3)/50*1000)
    print("Cathode voltage = {}".format(np.mean(s4)*500))  

    

    if do_spectrum:
      print("Max frequency, mode:")
      print("main mode = {}".format(m_main))
      print("main frequency = {}".format(max_f))   
    
    # add to lists
    run_list.append(run_num)
    I1.append(-np.mean(s1)/50*1000)
    I2.append(-np.mean(s2)/50*1000)
    I3.append(-np.mean(s3)/50*1000)
    I_tot_mean.append(np.abs(np.mean(Itot))*1000)
    I_tot.append(-Itot_main*1e3)
    V_cathode.append(Vc_main)
    V_cathode_mean.append(np.abs(np.mean(s4)))
    if do_spectrum:
      mode.append(m_main)
      main_frequency.append(f_main)
    if do_save:
      plt.close('all')




print("run #:")
print(run_list)
print("I1 mean:")
print(I1)
print("I2 mean:")
print(I2)
print("I3 mean:")
print(I3)
print("Itot main:")
print(I_tot)
print("Itot mean:")
print(I_tot_mean)
print("cathode voltage main:")
print(V_cathode)
print("cathode voltage mean:")
print(V_cathode_mean)
print("main mode:")
print(mode)
print("main frequency:")
print(main_frequency)

## Plots

I_main = np.array(I3)
I_dis = np.array(I_tot)
V_dis = np.array(V_cathode) - I_dis*65/1000

#plt.plot(I_dis, V_dis, 'ok')
#plt.show()



if do_show:
  plt.show()
