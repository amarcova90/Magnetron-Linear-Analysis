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
h_size_error = 5.5
v_size_error = 6

linesize = 3.0
linesize_large = 4.0

# subset of modes
unique_m = (3, 4)
do_decr = True
do_incr = False
do_show = True
do_save = False
do_spectrum = True

date = "2021-07-25"
date15 = "2021-07-23"
file_pramble = "WFM"
run_all = False
#runs_subset = [1]
runs_subset = range(116,138+1)
excelfile = "experiment/Testing_2020_new.xlsx"
voltage_tablename = "Plasma dis voltage"

max_i = 9.5

df = pd.read_excel(excelfile, sheet_name = date)

#### d = 2mm; p = 300mTorr

V_20_300 = []
I_20_300 = []
mode_20_300 = []
frequency_20_300 = []


run_l = 1
run_h = 26
# start routine
for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_20_300.append(df[voltage_tablename][row])
                    I_20_300.append(df["Plasma current (oscilloscope)"][row])
                    mode_20_300.append(df["mode"][row])
                    frequency_20_300.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_20_300.append(df[voltage_tablename][row])
                I_20_300.append(df["Plasma current (oscilloscope)"][row])
                mode_20_300.append(-100.0)
                frequency_20_300.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_20_300.append(df[voltage_tablename][row])
                I_20_300.append(df["Plasma current (oscilloscope)"][row])
                mode_20_300.append(100.0)
                frequency_20_300.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_20_300.append(df[voltage_tablename][row])
                I_20_300.append(df["Plasma current (oscilloscope)"][row])
                mode_20_300.append(0)
                frequency_20_300.append(-100.0)
            else:
                exit("ERROR 101")

V_20_300 = np.array(V_20_300)
I_20_300 = np.array(I_20_300)
mode_20_300 = np.array(mode_20_300)
frequency_20_300 = np.array(frequency_20_300)

##### d = 2mm; p = 200mTorr

V_20_200 = []
I_20_200 = []
mode_20_200 = []
frequency_20_200 = []

run_l = 27
run_h = 50
# start routine

for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_20_200.append(df[voltage_tablename][row])
                    I_20_200.append(df["Plasma current (oscilloscope)"][row])
                    mode_20_200.append(df["mode"][row])
                    frequency_20_200.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_20_200.append(df[voltage_tablename][row])
                I_20_200.append(df["Plasma current (oscilloscope)"][row])
                mode_20_200.append(-100.0)
                frequency_20_200.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_20_200.append(df[voltage_tablename][row])
                I_20_200.append(df["Plasma current (oscilloscope)"][row])
                mode_20_200.append(100.0)     
                frequency_20_200.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_20_200.append(df[voltage_tablename][row])
                I_20_200.append(df["Plasma current (oscilloscope)"][row])
                mode_20_200.append(0)
                frequency_20_200.append(-100.0)
            else:
                exit("ERROR 101")

V_20_200 = np.array(V_20_200)
I_20_200 = np.array(I_20_200)
mode_20_200 = np.array(mode_20_200)
frequency_20_200 = np.array(frequency_20_200)

##### d = 2mm; p = 100mTorr

V_20_100 = []
I_20_100 = []
mode_20_100 = []
frequency_20_100 = []

run_l = 51
run_h = 68
# start routine

for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_20_100.append(df[voltage_tablename][row])
                    I_20_100.append(df["Plasma current (oscilloscope)"][row])
                    mode_20_100.append(df["mode"][row])
                    frequency_20_100.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_20_100.append(df[voltage_tablename][row])
                I_20_100.append(df["Plasma current (oscilloscope)"][row])
                mode_20_100.append(-100.0)
                frequency_20_100.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_20_100.append(df[voltage_tablename][row])
                I_20_100.append(df["Plasma current (oscilloscope)"][row])
                mode_20_100.append(100.0)
                frequency_20_100.append(-100.0)
                         
            elif (df["mode"][row] == "mixed"):
                V_20_100.append(df[voltage_tablename][row])
                I_20_100.append(df["Plasma current (oscilloscope)"][row])
                mode_20_100.append(0)
                frequency_20_100.append(-100.0)
            else:
                exit("ERROR 101")

V_20_100 = np.array(V_20_100)
I_20_100 = np.array(I_20_100)
mode_20_100 = np.array(mode_20_100)
frequency_20_100 = np.array(frequency_20_100)

##### d = 2.5mm; p = 300mTorr

V_25_300 = []
I_25_300 = []
mode_25_300 = []
frequency_25_300 = []


run_l = 69
run_h = 96
# start routine

for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_25_300.append(df[voltage_tablename][row])
                    I_25_300.append(df["Plasma current (oscilloscope)"][row])
                    mode_25_300.append(df["mode"][row])
                    frequency_25_300.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_25_300.append(df[voltage_tablename][row])
                I_25_300.append(df["Plasma current (oscilloscope)"][row])
                mode_25_300.append(-100.0)
                frequency_25_300.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_25_300.append(df[voltage_tablename][row])
                I_25_300.append(df["Plasma current (oscilloscope)"][row])
                mode_25_300.append(100.0) 
                frequency_25_300.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_25_300.append(df[voltage_tablename][row])
                I_25_300.append(df["Plasma current (oscilloscope)"][row])
                mode_25_300.append(0)
                frequency_25_300.append(-100.0)
            else:
                exit("ERROR 101")
V_25_300 = np.array(V_25_300)
I_25_300 = np.array(I_25_300)
mode_25_300 = np.array(mode_25_300)
frequency_25_300 = np.array(frequency_25_300)


##### d = 2.5mm; p = 200mTorr

V_25_200 = []
I_25_200 = []
mode_25_200 = []
frequency_25_200 = []


run_l = 97
run_h = 115
# start routine

for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_25_200.append(df[voltage_tablename][row])
                    I_25_200.append(df["Plasma current (oscilloscope)"][row])
                    mode_25_200.append(df["mode"][row])
                    frequency_25_200.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_25_200.append(df[voltage_tablename][row])
                I_25_200.append(df["Plasma current (oscilloscope)"][row])
                mode_25_200.append(-100.0)
                frequency_25_200.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_25_200.append(df[voltage_tablename][row])
                I_25_200.append(df["Plasma current (oscilloscope)"][row])
                mode_25_200.append(100.0)
                frequency_25_200.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_25_200.append(df[voltage_tablename][row])
                I_25_200.append(df["Plasma current (oscilloscope)"][row])
                mode_25_200.append(0)
                frequency_25_200.append(-100.0)
            else:
                exit("ERROR 101")
V_25_200 = np.array(V_25_200)
I_25_200 = np.array(I_25_200)
mode_25_200 = np.array(mode_25_200)
frequency_25_200 = np.array(frequency_25_200)


##### d = 2.5mm; p = 100mTorr

V_25_100 = []
I_25_100 = []
mode_25_100 = []
frequency_25_100 = []

run_l = 116
run_h = 138
# start routine
for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_25_100.append(df[voltage_tablename][row])
                    I_25_100.append(df["Plasma current (oscilloscope)"][row])
                    mode_25_100.append(df["mode"][row])
                    frequency_25_100.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_25_100.append(df[voltage_tablename][row])
                I_25_100.append(df["Plasma current (oscilloscope)"][row])
                mode_25_100.append(-100.0)
                frequency_25_100.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_25_100.append(df[voltage_tablename][row])
                I_25_100.append(df["Plasma current (oscilloscope)"][row])
                mode_25_100.append(100.0)     
                frequency_25_100.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_25_100.append(df[voltage_tablename][row])
                I_25_100.append(df["Plasma current (oscilloscope)"][row])
                mode_25_100.append(0)
                frequency_25_100.append(-100.0)
            else:
                exit("ERROR 101")
V_25_100 = np.array(V_25_100)
I_25_100 = np.array(I_25_100)
mode_25_100 = np.array(mode_25_100)
frequency_25_100 = np.array(frequency_25_100)

## 1.5 mm

df = pd.read_excel(excelfile, sheet_name = date15)

##### d = 1.5mm; p = 300mTorr

V_15_300 = []
I_15_300 = []
mode_15_300 = []
frequency_15_300 = []


run_l = 39
run_h = 53
# start routine
for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_15_300.append(df[voltage_tablename][row])
                    I_15_300.append(df["Plasma current (oscilloscope)"][row])
                    mode_15_300.append(df["mode"][row])
                    frequency_15_300.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_15_300.append(df[voltage_tablename][row])
                I_15_300.append(df["Plasma current (oscilloscope)"][row])
                mode_15_300.append(-100.0)
                frequency_15_300.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_15_300.append(df[voltage_tablename][row])
                I_15_300.append(df["Plasma current (oscilloscope)"][row])
                mode_15_300.append(100.0)     
                frequency_15_300.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_15_300.append(df[voltage_tablename][row])
                I_15_300.append(df["Plasma current (oscilloscope)"][row])
                mode_15_300.append(0)
                frequency_15_300.append(-100.0)
            else:
                exit("ERROR 101")
V_15_300 = np.array(V_15_300)
I_15_300 = np.array(I_15_300)
mode_15_300 = np.array(mode_15_300)
frequency_15_300 = np.array(frequency_15_300)


##### d = 1.5mm; p = 200mTorr

V_15_200 = []
I_15_200 = []
mode_15_200 = []
frequency_15_200 = []


run_l = 54
run_h = 67
# start routine
for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :   
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_15_200.append(df[voltage_tablename][row])
                    I_15_200.append(df["Plasma current (oscilloscope)"][row])
                    mode_15_200.append(df["mode"][row])
                    frequency_15_200.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_15_200.append(df[voltage_tablename][row])
                I_15_200.append(df["Plasma current (oscilloscope)"][row])
                mode_15_200.append(-100.0)
                frequency_15_200.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_15_200.append(df[voltage_tablename][row])
                I_15_200.append(df["Plasma current (oscilloscope)"][row])
                mode_15_200.append(100.0)     
                frequency_15_200.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_15_200.append(df[voltage_tablename][row])
                I_15_200.append(df["Plasma current (oscilloscope)"][row])
                mode_15_200.append(0)
                frequency_15_200.append(-100.0)
            else:
                exit("ERROR 101")
V_15_200 = np.array(V_15_200)
I_15_200 = np.array(I_15_200)
mode_15_200 = np.array(mode_15_200)
frequency_15_200 = np.array(frequency_15_200)


##### d = 1.5mm; p = 100mTorr

V_15_100 = []
I_15_100 = []
mode_15_100 = []
frequency_15_100 = []


run_l = 68
run_h = 91
# start routine
for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run >= run_l and run <= run_h and df["Plasma current (oscilloscope)"][row] < max_i :
            print(run)
            if type(df["mode"][row]) is float:
                if np.isnan(df["mode"][row]) == False and df["plot"][row]:
                    V_15_100.append(df[voltage_tablename][row])
                    I_15_100.append(df["Plasma current (oscilloscope)"][row])
                    mode_15_100.append(df["mode"][row])
                    frequency_15_100.append(df["main frequency"][row])
                else:
                    continue
            elif (df["mode"][row] == "neg"):
                V_15_100.append(df[voltage_tablename][row])
                I_15_100.append(df["Plasma current (oscilloscope)"][row])
                mode_15_100.append(-100.0)
                frequency_15_100.append(-100.0)
            elif (df["mode"][row] == "pos"):
                V_15_100.append(df[voltage_tablename][row])
                I_15_100.append(df["Plasma current (oscilloscope)"][row])
                mode_15_100.append(100.0)
                frequency_15_100.append(-100.0)
            elif (df["mode"][row] == "mixed"):
                V_15_100.append(df[voltage_tablename][row])
                I_15_100.append(df["Plasma current (oscilloscope)"][row])
                mode_15_100.append(0)
                frequency_15_100.append(-100.0)
            else:
                exit("ERROR 101")
V_15_100 = np.array(V_15_100)
I_15_100 = np.array(I_15_100)
mode_15_100 = np.array(mode_15_100)
frequency_15_100 = np.array(frequency_15_100)


# treshold current error bars

It_15_300 = np.array([4.0239, 5.366])
It_15_200 = np.array([5.8205, 6.7011])
It_15_100 = np.array([])
It_20_300 = np.array([3.0407, 3.4876])
It_20_200 = np.array([4.9101, 6.0614])
It_20_100 = np.array([])
It_25_300 = np.array([2.2865, 2.4535])
It_25_200 = np.array([3.4879, 4.487])
It_25_100 = np.array([5.0823, 5.5601])

It_mean_15_300 = np.mean(It_15_300)
It_error_15_300 = np.abs(It_15_300-It_mean_15_300)
It_mean_20_300 = np.mean(It_20_300)
It_error_20_300 = np.abs(It_20_300-It_mean_20_300)
It_mean_25_300 = np.mean(It_25_300)
It_error_25_300 = np.abs(It_25_300-It_mean_25_300)

It_mean_300_v = np.array([It_mean_15_300, It_mean_20_300, It_mean_25_300])
It_error_300_v = np.transpose(np.array([It_error_15_300, It_error_20_300, It_error_25_300]))

It_mean_15_200 = np.mean(It_15_200)
It_error_15_200 = np.abs(It_15_200-It_mean_15_200)
It_mean_20_200 = np.mean(It_20_200)
It_error_20_200 = np.abs(It_20_200-It_mean_20_200)
It_mean_25_200 = np.mean(It_25_200)
It_error_25_200 = np.abs(It_25_200-It_mean_25_200)

It_mean_200_v = np.array([It_mean_15_200, It_mean_20_200, It_mean_25_200])
It_error_200_v = np.transpose(np.array([It_error_15_200, It_error_20_200, It_error_25_200]))

It_mean_25_100= np.mean(It_25_100)
It_error_25_100 = np.abs(It_25_100-It_mean_25_100)

It_mean_100_v = np.array([It_mean_25_100])
It_error_100_v = np.transpose(np.array([It_error_25_100]))


d_v = np.array([1.5, 2., 2.5])


plt.figure(figsize=(h_size_error, v_size_error))
plt.errorbar(d_v, It_mean_300_v, yerr = It_error_300_v, fmt ='.k', markersize = 0, capsize = 20)
#plt.xlabel("distance (mm)")
#plt.ylabel("treshold current (mA)")
#plt.ylim((0.0, 8.0))
#plt.tight_layout()
#plt.savefig("figures_paper2/errorbars_p=300.png")


#plt.figure(figsize=(h_size*0.8, v_size*0.8))
plt.errorbar(d_v, It_mean_200_v, yerr = It_error_200_v, fmt ='.b', markersize = 0, capsize = 20)

plt.errorbar(np.array([2.5]), It_mean_100_v, yerr = It_error_100_v, fmt ='.r', markersize = 0, capsize = 20)

plt.xlabel("distance (mm)")
plt.ylabel("threshold current (mA)")
plt.legend(("300mTorr","200mTorr", "100mTorr"))
plt.xlim([1.2, 2.8])
plt.ylim((0.0, 8.0))
plt.tight_layout()
plt.savefig("figures_paper2/errorbars.png")


######################## mode theory figures #################

####  m = 6,7; d = 2.5mm; p = 200 mTorr: 2021-08-18


date = "2021-08-18"
excelfile = "experiment/Testing_2020_new3.xlsx"
df = pd.read_excel(excelfile, sheet_name = date)

# m = 6
runs_sub = (12, 15, 17, 23)

V_6 = []
I_6 = []
f_6 = []


for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run in runs_sub:
            V_6.append(df["Plasma dis voltage"][row])
            I_6.append(df["Plasma current (oscilloscope)"][row])
            f_6.append(df["main frequency"][row])

            
V_6 = np.array(V_6)
I_6 = np.array(I_6)
f_6 = np.array(f_6)

# m = 7
runs_sub = (10, 18, 21)

V_7 = []
I_7 = []
f_7 = []


for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run in runs_sub:
            V_7.append(df["Plasma dis voltage"][row])
            I_7.append(df["Plasma current (oscilloscope)"][row])
            f_7.append(df["main frequency"][row])

            
V_7 = np.array(V_7)
I_7 = np.array(I_7)
f_7 = np.array(f_7)




####  m = 5; d = 2 mm; p = 200 mTorr: 2021-09-15


date = "2021-09-15"
excelfile = "experiment/Testing_2020_new3.xlsx"
df = pd.read_excel(excelfile, sheet_name = date)

# m = 6
runs_sub = (2, 3, 4, 5, 6)

V_5 = []
I_5 = []
f_5 = []


for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run in runs_sub:
            V_5.append(df["Plasma dis voltage"][row])
            I_5.append(df["Plasma current (oscilloscope)"][row])
            f_5.append(df["main frequency"][row])

            
V_5 = np.array(V_5)
I_5 = np.array(I_5)
f_5 = np.array(f_5)

# m = 4
runs_sub = (17, 19, 20, 21, 23)

V_4 = []
I_4 = []
f_4 = []


for row,run in enumerate(df["Run"]):
    if (type(run) is int):
        if run in runs_sub:
            V_4.append(df["Plasma dis voltage"][row])
            I_4.append(df["Plasma current (oscilloscope)"][row])
            f_4.append(df["main frequency"][row])

            
V_4 = np.array(V_4)
I_4 = np.array(I_4)
f_4 = np.array(f_4)


plt.figure(figsize=(h_size, v_size))
plt.plot(V_4, f_4, "^b")
plt.plot(V_5, f_5, "or")
plt.legend(("m = 4", "m = 5"))
plt.xlabel("discharge voltage (V)")
plt.ylabel("frequency (MHz)")
plt.title("d = 2mm")
plt.tight_layout()
plt.savefig("figures_paper2/Vf_d=2.png")

plt.figure(figsize=(h_size, v_size))
plt.plot(I_4, f_4, "^b")
plt.plot(I_5, f_5, "or")
plt.legend(("m = 4", "m = 5"))
plt.xlabel("discharge current (mA)")
plt.ylabel("frequency (MHz)")
plt.title("d = 2mm")
plt.tight_layout()
plt.savefig("figures_paper2/if_d=2.png")

# plt.figure(figsize=(h_size, v_size))
# plt.plot(V_6, f_6, "sk")
# plt.plot(V_7, f_7, "Dm")
# plt.legend(("m = 6", "m = 7"))
# # plt.ylim([2.25,3])
# # plt.xlim([268,271])
# plt.xlabel("discharge voltage (V)")
# plt.ylabel("frequency (MHz)")
# plt.title("d = 2.5 mm")
# plt.tight_layout()
# plt.savefig("figures_paper2/Vf_d=25.png")

# plt.figure(figsize=(h_size, v_size))
# plt.plot(I_6, f_6, "sk")
# plt.plot(I_7, f_7, "Dm")
# plt.legend(("m = 4", "m = 5"))
# plt.xlabel("discharge current (I)")
# plt.ylabel("frequency (MHz)")
# plt.title("d = 2mm")
# plt.tight_layout()

# plt.figure(figsize=(h_size, v_size))
# plt.plot(V_4, f_4, "^b")
# plt.plot(V_5, f_5, "or")
# plt.plot(V_6, f_6, "sk")
# plt.plot(V_7, f_7, "Dm")
# plt.legend(("m = 4", "m = 5", "m = 6", "m = 7"))
# plt.xlabel("discharge voltage (V)")
# plt.ylabel("frequency (MHz)")
# plt.tight_layout()
# plt.savefig("figures_paper2/Vf_both.png")


##################### IV and VI figures #############################

nExB_color = "r"
nExB_marker = "<"

ExB_color = "g"
ExB_marker = ">"

mixed_color = "b"
mixed_marker = "o"

black_line1 = "--k"
black_line2 = "-.k"
black_line3 = ":k"

ylimits = [245, 300]
xlimits = [1.5, 10.0]

marker_size = 8


## IV 20 

plt.figure(figsize=(h_size, v_size))

plt.plot(V_20_300 ,I_20_300, black_line1, linewidth = 1, label='_nolegend_')
plt.plot(V_20_300[mode_20_300 > 0.0] ,I_20_300[mode_20_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(V_20_300[mode_20_300 < 0.0] ,I_20_300[mode_20_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(V_20_300[mode_20_300 == 0.0] ,I_20_300[mode_20_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_20_200 ,I_20_200, black_line2, linewidth = 1, label='_nolegend_')
plt.plot(V_20_200[mode_20_200 > 0.0] ,I_20_200[mode_20_200 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_20_200[mode_20_200 < 0.0] ,I_20_200[mode_20_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_20_200[mode_20_200 == 0.0] ,I_20_200[mode_20_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_20_100 ,I_20_100, black_line3, linewidth = 1, label='_nolegend_')
plt.plot(V_20_100[mode_20_100 > 0.0] ,I_20_100[mode_20_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_20_100[mode_20_100 < 0.0] ,I_20_100[mode_20_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_20_100[mode_20_100 == 0.0] ,I_20_100[mode_20_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.xlabel("discharge voltage (V)")
plt.ylabel("discharge current (mA)")
plt.title("d = 2mm")
plt.ylim(xlimits)
plt.xlim(ylimits)
plt.tight_layout()
plt.savefig("figures_paper2/IV_d=20.png")

## VI 20 

plt.figure(figsize=(v_size, h_size))

plt.plot(I_20_300, V_20_300 , black_line1, linewidth = 1, label='_nolegend_')
plt.plot(I_20_300[mode_20_300 > 0.0], V_20_300[mode_20_300 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(I_20_300[mode_20_300 < 0.0], V_20_300[mode_20_300 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(I_20_300[mode_20_300 == 0.0], V_20_300[mode_20_300 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_20_200, V_20_200 , black_line2, linewidth = 1, label='_nolegend_')
plt.plot(I_20_200[mode_20_200 > 0.0], V_20_200[mode_20_200 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_20_200[mode_20_200 < 0.0], V_20_200[mode_20_200 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_20_200[mode_20_200 == 0.0], V_20_200[mode_20_200 == 0.0] ,mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_20_100, V_20_100, black_line3, linewidth = 1, label='_nolegend_')
plt.plot(I_20_100[mode_20_100 > 0.0], V_20_100[mode_20_100 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_20_100[mode_20_100 < 0.0], V_20_100[mode_20_100 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_20_100[mode_20_100 == 0.0], V_20_100[mode_20_100 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color)
#plt.legend(("- ExB","  ExB", "mixed"))

plt.ylabel("discharge voltage (V)")
plt.xlabel("discharge current (mA)")
plt.title("d = 2mm")
plt.ylim(ylimits)
plt.xlim(xlimits)
plt.tight_layout()
plt.savefig("figures_paper2/VI_d=20.png")

## IV 25

plt.figure(figsize=(h_size, v_size))

plt.plot(V_25_300 ,I_25_300, black_line1, linewidth = 1, label='_nolegend_')
plt.plot(V_25_300[mode_25_300 > 0.0] ,I_25_300[mode_25_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(V_25_300[mode_25_300 < 0.0] ,I_25_300[mode_25_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(V_25_300[mode_25_300 == 0.0] ,I_25_300[mode_25_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_25_200 ,I_25_200, black_line2, linewidth = 1, label='_nolegend_')
plt.plot(V_25_200[mode_25_200 > 0.0] ,I_25_200[mode_25_200 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_25_200[mode_25_200 < 0.0] ,I_25_200[mode_25_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_25_200[mode_25_200 == 0.0] ,I_25_200[mode_25_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_25_100 ,I_25_100, black_line3, linewidth = 1, label='_nolegend_')
plt.plot(V_25_100[mode_25_100 > 0.0] ,I_25_100[mode_25_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_25_100[mode_25_100 < 0.0] ,I_25_100[mode_25_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_25_100[mode_25_100 == 0.0] ,I_25_100[mode_25_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.legend(("- ExB","  ExB", "mixed"))

plt.xlabel("discharge voltage (V)")
plt.ylabel("discharge current (mA)")
plt.title("d = 2.5mm")
plt.ylim(xlimits)
plt.xlim(ylimits)
plt.tight_layout()
plt.savefig("figures_paper2/IV_d=25.png")


## VI 25

plt.figure(figsize=(v_size, h_size))

plt.plot(I_25_300, V_25_300 , black_line1, linewidth = 1, label='_nolegend_')
plt.plot(I_25_300[mode_25_300 > 0.0], V_25_300[mode_25_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(I_25_300[mode_25_300 < 0.0], V_25_300[mode_25_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(I_25_300[mode_25_300 == 0.0], V_25_300[mode_25_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_25_200, V_25_200 , black_line2, linewidth = 1, label='_nolegend_')
plt.plot(I_25_200[mode_25_200 > 0.0], V_25_200[mode_25_200 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_25_200[mode_25_200 < 0.0], V_25_200[mode_25_200 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_25_200[mode_25_200 == 0.0], V_25_200[mode_25_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_25_100, V_25_100 , black_line3, linewidth = 1, label='_nolegend_')
plt.plot(I_25_100[mode_25_100 > 0.0], V_25_100[mode_25_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_25_100[mode_25_100 < 0.0], V_25_100[mode_25_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_25_100[mode_25_100 == 0.0], V_25_100[mode_25_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.legend(("- ExB","  ExB", "mixed"))

plt.ylabel("discharge voltage (V)")
plt.xlabel("discharge current (mA)")
plt.title("d = 2.5mm")
plt.ylim(ylimits)
plt.xlim(xlimits)
plt.tight_layout()
plt.savefig("figures_paper2/VI_d=25.png")

## IV 15

plt.figure(figsize=(h_size, v_size))

plt.plot(V_15_300 ,I_15_300, black_line1, linewidth = 1, label='_nolegend_')
plt.plot(V_15_300[mode_15_300 > 0.0] ,I_15_300[mode_15_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(V_15_300[mode_15_300 < 0.0] ,I_15_300[mode_15_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(V_15_300[mode_15_300 == 0.0] ,I_15_300[mode_15_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_15_200[np.argsort(I_15_200)], I_15_200[np.argsort(I_15_200)], black_line2, linewidth = 1, label='_nolegend_')
plt.plot(V_15_200[mode_15_200 > 0.0] ,I_15_200[mode_15_200 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_15_200[mode_15_200 < 0.0] ,I_15_200[mode_15_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_15_200[mode_15_200 == 0.0] ,I_15_200[mode_15_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_15_100[np.argsort(I_15_100)], I_15_100[np.argsort(I_15_100)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(V_15_100[mode_15_100 > 0.0] ,I_15_100[mode_15_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_15_100[mode_15_100 < 0.0] ,I_15_100[mode_15_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_15_100[mode_15_100 == 0.0] ,I_15_100[mode_15_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

#plt.legend(("- ExB","  ExB", "mixed"))

plt.xlabel("discharge voltage (V)")
plt.ylabel("discharge current (mA)")
plt.title("d = 1.5mm")
plt.ylim(xlimits)
plt.xlim(ylimits)
plt.tight_layout()
plt.savefig("figures_paper2/IV_d=15.png")


## VI 15

plt.figure(figsize=(v_size, h_size))

plt.plot(I_15_300, V_15_300, black_line1, linewidth = 1, label='_nolegend_')
plt.plot(I_15_300[mode_15_300 > 0.0], V_15_300[mode_15_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(I_15_300[mode_15_300 < 0.0], V_15_300[mode_15_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(I_15_300[mode_15_300 == 0.0], V_15_300[mode_15_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_15_200[np.argsort(I_15_200)], V_15_200[np.argsort(I_15_200)], black_line2, linewidth = 1, label='_nolegend_')
plt.plot(I_15_200[mode_15_200 > 0.0], V_15_200[mode_15_200 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_15_200[mode_15_200 < 0.0], V_15_200[mode_15_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_15_200[mode_15_200 == 0.0], V_15_200[mode_15_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_15_100[np.argsort(I_15_100)], V_15_100[np.argsort(I_15_100)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(I_15_100[mode_15_100 > 0.0], V_15_100[mode_15_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_15_100[mode_15_100 < 0.0], V_15_100[mode_15_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_15_100[mode_15_100 == 0.0], V_15_100[mode_15_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

#plt.legend(("- ExB","  ExB", "mixed"))

plt.ylabel("discharge voltage (V)")
plt.xlabel("discharge current (mA)")
plt.title("d = 1.5mm")
plt.xlim(xlimits)
plt.ylim(ylimits)
plt.tight_layout()
plt.savefig("figures_paper2/VI_d=15.png")


## same pressure different distances

# IV 100

plt.figure(figsize=(h_size, v_size))

plt.plot(V_15_100[np.argsort(I_15_100)], I_15_100[np.argsort(I_15_100)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(V_15_100[mode_15_100 > 0.0] ,I_15_100[mode_15_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(V_15_100[mode_15_100 < 0.0] ,I_15_100[mode_15_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(V_15_100[mode_15_100 == 0.0] ,I_15_100[mode_15_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_20_100 ,I_20_100, black_line2, linewidth = 1, label='_nolegend_')
plt.plot(V_20_100[mode_20_100 > 0.0] ,I_20_100[mode_20_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_20_100[mode_20_100 < 0.0] ,I_20_100[mode_20_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_20_100[mode_20_100 == 0.0] ,I_20_100[mode_20_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')

plt.plot(V_25_100 ,I_25_100, black_line3, linewidth = 1, label='_nolegend_')
plt.plot(V_25_100[mode_25_100 > 0.0] ,I_25_100[mode_25_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_25_100[mode_25_100 < 0.0] ,I_25_100[mode_25_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_25_100[mode_25_100 == 0.0] ,I_25_100[mode_25_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')


#plt.legend(("- ExB","  ExB", "mixed"))

plt.xlabel("discharge voltage (V)")
plt.ylabel("discharge current (mA)")
plt.title("p = 100mTorr")
plt.tight_layout()
plt.savefig("figures_paper2/IV_p=100.png")


# VI 100

plt.figure(figsize=(v_size, h_size))

plt.plot(I_15_100[np.argsort(I_15_100)], V_15_100[np.argsort(I_15_100)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(I_15_100[mode_15_100 > 0.0], V_15_100[mode_15_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(I_15_100[mode_15_100 < 0.0], V_15_100[mode_15_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(I_15_100[mode_15_100 == 0.0], V_15_100[mode_15_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_20_100, V_20_100 , black_line2, linewidth = 1, label='_nolegend_')
plt.plot(I_20_100[mode_20_100 > 0.0], V_20_100[mode_20_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_20_100[mode_20_100 < 0.0], V_20_100[mode_20_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_20_100[mode_20_100 == 0.0], V_20_100[mode_20_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')

plt.plot(I_25_100, V_25_100, black_line3, linewidth = 1, label='_nolegend_')
plt.plot(I_25_100[mode_25_100 > 0.0], V_25_100[mode_25_100 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_25_100[mode_25_100 < 0.0], V_25_100[mode_25_100 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_25_100[mode_25_100 == 0.0], V_25_100[mode_25_100 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')


#plt.legend(("- ExB","  ExB", "mixed"))

plt.ylabel("discharge voltage (V)")
plt.xlabel("discharge current (mA)")
plt.title("p = 100mTorr")
plt.tight_layout()
plt.savefig("figures_paper2/VI_p=100.png")

## IV 200

plt.figure(figsize=(h_size, v_size))

plt.plot(V_15_200[np.argsort(I_15_200)], I_15_200[np.argsort(I_15_200)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(V_15_200[mode_15_200 > 0.0] ,I_15_200[mode_15_200 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(V_15_200[mode_15_200 < 0.0] ,I_15_200[mode_15_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(V_15_200[mode_15_200 == 0.0] ,I_15_200[mode_15_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_20_200 ,I_20_200, black_line2, linewidth = 1, label='_nolegend_')
plt.plot(V_20_200[mode_20_200 > 0.0] ,I_20_200[mode_20_200 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_20_200[mode_20_200 < 0.0] ,I_20_200[mode_20_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_20_200[mode_20_200 == 0.0] ,I_20_200[mode_20_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')

plt.plot(V_25_200 ,I_25_200, black_line3, linewidth = 1, label='_nolegend_')
plt.plot(V_25_200[mode_25_200 > 0.0] ,I_25_200[mode_25_200 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_25_200[mode_25_200 < 0.0] ,I_25_200[mode_25_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_25_200[mode_25_200 == 0.0] ,I_25_200[mode_25_200 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')


#plt.legend(("- ExB","  ExB", "mixed"))

plt.xlabel("discharge voltage (V)")
plt.ylabel("discharge current (mA)")
plt.title("p = 200mTorr")
plt.tight_layout()
plt.savefig("figures_paper2/IV_p=200.png")

# VI 200

plt.figure(figsize=(v_size, h_size))

plt.plot(I_15_200[np.argsort(I_15_200)], V_15_200[np.argsort(I_15_200)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(I_15_200[mode_15_200 > 0.0],V_15_200[mode_15_200 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(I_15_200[mode_15_200 < 0.0], V_15_200[mode_15_200 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(I_15_200[mode_15_200 == 0.0], V_15_200[mode_15_200 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_20_200, V_20_200 , black_line2, linewidth = 1, label='_nolegend_')
plt.plot(I_20_200[mode_20_200 > 0.0], V_20_200[mode_20_200 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_20_200[mode_20_200 < 0.0], V_20_200[mode_20_200 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_20_200[mode_20_200 == 0.0], V_20_200[mode_20_200 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')

plt.plot(I_25_200, V_25_200 , black_line3, linewidth = 1, label='_nolegend_')
plt.plot(I_25_200[mode_25_200 > 0.0], V_25_200[mode_25_200 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_25_200[mode_25_200 < 0.0], V_25_200[mode_25_200 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_25_200[mode_25_200 == 0.0], V_25_200[mode_25_200 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')


#plt.legend(("- ExB","  ExB", "mixed"))

plt.ylabel("discharge voltage (V)")
plt.xlabel("discharge current (mA)")
plt.title("p = 200mTorr")
plt.tight_layout()
plt.savefig("figures_paper2/VI_p=200.png")


# IV 300

plt.figure(figsize=(h_size, v_size))

plt.plot(V_15_300[np.argsort(I_15_300)], I_15_300[np.argsort(I_15_300)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(V_15_300[mode_15_300 > 0.0] ,I_15_300[mode_15_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(V_15_300[mode_15_300 < 0.0] ,I_15_300[mode_15_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(V_15_300[mode_15_300 == 0.0] ,I_15_300[mode_15_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(V_20_300 ,I_20_300, black_line2, linewidth = 1, label='_nolegend_')
plt.plot(V_20_300[mode_20_300 > 0.0] ,I_20_300[mode_20_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_20_300[mode_20_300 < 0.0] ,I_20_300[mode_20_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_20_300[mode_20_300 == 0.0] ,I_20_300[mode_20_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')

plt.plot(V_25_300 ,I_25_300, black_line3, linewidth = 1, label='_nolegend_')
plt.plot(V_25_300[mode_25_300 > 0.0] ,I_25_300[mode_25_300 > 0.0], nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(V_25_300[mode_25_300 < 0.0] ,I_25_300[mode_25_300 < 0.0], ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(V_25_300[mode_25_300 == 0.0] ,I_25_300[mode_25_300 == 0.0], mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')


plt.legend(("- ExB","  ExB", "mixed"))

plt.xlabel("discharge voltage (V)")
plt.ylabel("discharge current (mA)")
plt.title("p = 300mTorr")
plt.tight_layout()
plt.savefig("figures_paper2/IV_p=300.png")


# VI 300

plt.figure(figsize=(v_size, h_size))

plt.plot(I_15_300[np.argsort(I_15_300)], V_15_300[np.argsort(I_15_300)], black_line1, linewidth = 1, label='_nolegend_')
plt.plot(I_15_300[mode_15_300 > 0.0], V_15_300[mode_15_300 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color)
plt.plot(I_15_300[mode_15_300 < 0.0], V_15_300[mode_15_300 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color)
plt.plot(I_15_300[mode_15_300 == 0.0], V_15_300[mode_15_300 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color)

plt.plot(I_20_300, V_20_300 , black_line2, linewidth = 1, label='_nolegend_')
plt.plot(I_20_300[mode_20_300 > 0.0], V_20_300[mode_20_300 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_20_300[mode_20_300 < 0.0], V_20_300[mode_20_300 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_20_300[mode_20_300 == 0.0], V_20_300[mode_20_300 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')

plt.plot(I_25_300, V_25_300 , black_line3, linewidth = 1, label='_nolegend_')
plt.plot(I_25_300[mode_25_300 > 0.0], V_25_300[mode_25_300 > 0.0] , nExB_marker, markersize = marker_size, color = nExB_color, label='_nolegend_')
plt.plot(I_25_300[mode_25_300 < 0.0], V_25_300[mode_25_300 < 0.0] , ExB_marker, markersize = marker_size, color = ExB_color, label='_nolegend_')
plt.plot(I_25_300[mode_25_300 == 0.0], V_25_300[mode_25_300 == 0.0] , mixed_marker, markersize = marker_size, color = mixed_color, label='_nolegend_')


plt.legend(("- ExB","  ExB", "mixed"))

plt.ylabel("discharge voltage (V)")
plt.xlabel("discharge current (mA)")
plt.title("p = 300mTorr")
plt.tight_layout()
plt.savefig("figures_paper2/VI_p=300.png")


######### wave dispertion




if do_show:
    plt.show()
