from scipy import signal
import numpy as np
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

# s1_file   = "experiment/Ar/07252020/WF_P13.csv"
# s2_file   = "experiment/Ar/07252020/WF_P14.csv"
# s3_file   = "experiment/Ar/07252020/WF_P15.csv"
# s4_file   = "experiment/Ar/07252020/WF_P16.csv"

s1_file   = "experiment/Ar/01102021/WF_P07.csv"
s2_file   = "experiment/Ar/01102021/WF_P09.csv"
s3_file   = "experiment/Ar/01102021/WF_P09.csv"
s4_file   = "experiment/Ar/01102021/WF_P10.csv"

# s1_file   = "experiment/Ar/WF_P04.csv"
# s2_file   = "experiment/Ar/WF_P05.csv"
# s3_file   = "experiment/Ar/WF_P06.csv"
# s4_file   = "experiment/Ar/07252020/WF_P16.csv"

gas = s1_file[11:13]
date = s1_file[14:22]
s1_name = s1_file[23:-4]

# make directory to save figures
cwd = os.getcwd()

# save the figures in the data folder
# path = os.path.join(cwd,"experiment",gas,date,"figures") 

# # create directory if it doesn't exist
# if not os.path.isdir(path):
  # os.mkdir(path)



#df = pd.read_excel(os.path.join(cwd,"experiment","Testing.xlsx"), sheet_name = "July 25 2020") 

#print(df.index)

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

R0 = 2.5e-3
#dx = 10.9*np.pi/180.0
dx = 10*np.pi/180.0
dx = 20*np.pi/180.0
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

print("main mode = {}".format(kx12[indx]))
print("main frequency = {}".format(f_v[indx]))

plt.figure(figsize=(h_size, v_size))
plt.plot(f_v,kx12)
  

## spectrum analysis

nf = 1000
nk = 1000
n = len(s1_detrended)

kNyq = 2*np.pi/dx/2.0			# Nyquist wavenumber
fNyq = 1.0/dt/2.0		# Nyquist frequency
fupper = fNyq/2.001	# maximum resolvable frequency
flower = 4.0/dt/n		# minimum resolvable frequency

f = np.logspace(np.log10(flower),np.log10(fupper),nf)
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
    
    
print("ballast voltage drop = {}".format(np.mean(s4)*500))     

print(S)

plt.figure(figsize=(h_size, v_size))
plt.pcolormesh(k,fz*1e-3,np.log10(S), cmap = 'jet')
# plt.pcolormesh(k + 2*kNyq,fz*1e-3,np.log10(S), cmap = 'jet')
# plt.pcolormesh(k - 2*kNyq,fz*1e-3,np.log10(S), cmap = 'jet')
plt.xlabel("mode")
plt.ylabel("frequency [kHz]")
plt.ylim((0,10*1e3))
plt.colorbar()
plt.tight_layout()
# filename = os.path.join(path,"pspect_pcolormesh.png") 
# plt.savefig(filename)  

# Something this looks better:
# plt.figure(figsize=(h_size, v_size))
# plt.contourf(k,fz*1e-3,np.log10(S), cmap = 'jet')
# plt.xlabel("mode")
# plt.ylabel("frequency [kHz]")
# plt.ylim((0,10*1e3))
# plt.colorbar()
# plt.tight_layout()
 
# filename = os.path.join(path,"pspect_contourf.png") 
# plt.savefig(filename)  

# print("f = {}".format(xfft))
#for i in range(len(ff1))


# fupper = fNyq/2.0
# flower = 4.0/dt/len(s1)
# f = np.linspace(flower,fupper,len(s1)/2.0)

plt.figure(figsize=(h_size, v_size))
plt.plot(f_v*1e-3,np.abs(ff1))
plt.plot(f_v*1e-3,np.abs(ff2))
plt.xlabel("frequency [kHz]")
plt.ylabel("signal amp [mA]")
plt.legend(["s1","s2"])
plt.tight_layout()
plt.xlim((0,10*1e3))
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
# filename = os.path.join(path,"FFT_logy.png") 
# plt.savefig(filename)  

plt.figure(figsize=(h_size, v_size))
plt.plot(t*1e6,s1_detrended)
plt.plot(t*1e6,s2_detrended)
plt.plot(t*1e6,s3_detrended)
plt.xlabel("time [$\mu$s]")
plt.ylabel("signal")
plt.legend(["s1","s2","s3"])
#plt.xlim((0,10))
plt.title("detrended signals")


plt.figure(figsize=(h_size, v_size))
plt.plot(t*1e6,s1)
plt.plot(t*1e6,s2)
plt.plot(t*1e6,s3)
plt.xlabel("time [$\mu$s]")
plt.ylabel("signal")
plt.legend(["s1","s2", "s3"])
#plt.xlim((0,10))
plt.title("raw signals")

# filename = os.path.join(path,"signals.png") 
# plt.savefig(filename) 

# plt.figure(figsize=(h_size, v_size))
# plt.plot(t*1e6,s1)
# plt.plot(t*1e6,s2)
# plt.plot(t*1e6,s3)
# #plt.plot(t*1e6,(s1+s2+s3)*1000)
# plt.xlabel("time [$\mu$s]")
# plt.ylabel("current [mA]")
# plt.legend(["s1","s2","s3"])

# plt.figure(figsize=(h_size, v_size))
# plt.plot(t*1e6,s4)
# #plt.plot(t*1e6,(s1+s2+s3)*1000)
# plt.xlabel("time [$\mu$s]")
# plt.ylabel("current [mA]")
# plt.legend(["s4"])



nf = len(ff1)

f = np.linspace(flower,fupper,nf) # array of frequencies
scale = 1.0/(f*dt) # array of wavelet scales


# for i in range(nf):
  # wx, xfft = pywt.cwt(s1,scale[i],'morl')
  # Z1[i] = wx

# plt.figure(figsize=(h_size, v_size))
# plt.contour(t,fz,np.real(Z1))

# plt.figure(figsize=(h_size, v_size))
# plt.plot(t*1e6,wx[0])

Z1 = signal.cwt(s1_detrended, signal.morlet2, scale) 

Z2 = signal.cwt(s2_detrended, signal.morlet2, scale) 

fz = pywt.scale2frequency('morl', scale)/dt

plt.figure(figsize=(h_size, v_size))
plt.contourf(t*1e6,np.log10(fz),np.abs(Z1), cmap = 'jet')
plt.colorbar()
plt.ylim((0,10*1e3))
plt.xlabel("time [$\mu$s]")
plt.ylabel("frequency [kHz]")
plt.title("signal 1")
plt.tight_layout()
# filename = os.path.join(path,"signal1Morlet_contourf.png") 
# plt.savefig(filename) 

plt.figure(figsize=(h_size, v_size))
plt.contourf(t*1e6,np.log10(fz),np.abs(Z2), cmap = 'jet')
plt.colorbar()
plt.ylim((0,10*1e3))
plt.xlabel("time [$\mu$s]")
plt.ylabel("frequency [kHz]")
plt.title("signal 1")
plt.tight_layout()
# filename = os.path.join(path,"signal2Morlet_contourf.png") 
# plt.savefig(filename) 


plt.figure(figsize=(h_size, v_size))
plt.contour(t,fz*1e-3,np.abs(Z1))
plt.colorbar()
plt.ylim((0,10*1e3))
plt.xlabel("time [$\mu$s]")
plt.ylabel("frequency [kHz]")
plt.title("signal 1")
plt.tight_layout()
# filename = os.path.join(path,"signal1Morlet_contour.png") 
# plt.savefig(filename) 


plt.figure(figsize=(h_size, v_size))
plt.contour(t,fz*1e-3,np.abs(Z2))
plt.colorbar()
plt.ylim((0,10*1e3))
plt.xlabel("time [$\mu$s]")
plt.ylabel("frequency [kHz]")
plt.title("signal 2")
plt.tight_layout() 
# filename = os.path.join(path,"signal2Morlet_contour.png") 
# plt.savefig(filename) 


plt.show()

