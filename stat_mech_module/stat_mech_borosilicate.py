# -*- coding: utf-8 -*-
"""
Created on Tue May 29 10:24:49 2018

@author: msb
"""
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import csv
import scipy.optimize
import math
import warnings
import os

warnings.filterwarnings("ignore", category=RuntimeWarning)


fil = 'Qn2'
fil2 = 'Bn3'

fil3 = "Qn_k2"
fil4 = "Bn_k2"

fil5 = "Qn_k4"
fil6 = "BnMD_k4"

fil7 = "Qn_k6"
fil8 = "K_BnMD_k6"

mod_si_data = []
mod_b_data = []

k_si_data = []
k_b_data = []
#N3_data = []
N4_data = []
#N2_data = []
#N1_data = []
#N0_data = []

Q4_data = []
Q3_data = []
Q2_data = []
Q1_data = []
Q0_data = []

if not os.path.exists(f"{fil}_data"):
        os.mkdir(f"{fil}_data")


with open(f"{fil}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        k_si_data.append(row[0])
        mod_si_data.append(row[1])
        Q4_data.append(row[2])
        Q3_data.append(row[3])
        Q2_data.append(row[4])
        Q1_data.append(row[5])
        Q0_data.append(row[6])

with open(f"{fil2}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        k_b_data.append(row[0])
        mod_b_data.append(row[1])
#        N3_data.append(row[2])
        N4_data.append(row[2])
#        N2_data.append(row[4])
#        N1_data.append(row[5])
#        N0_data.append(row[6])


mod_si_data = [float(i) for i in mod_si_data]
mod_b_data = [float(i) for i in mod_b_data]

k_si_data = [float(i) for i in k_si_data]
k_si_nr = list(range(len(k_si_data)))
k_b_data = [float(i)for i in k_b_data]
k_b_nr = list(range(len(k_b_data)))

#N3_data = [float(i) for i in N3_data]
N4_data = [float(i) for i in N4_data]
#N2_data = [float(i) for i in N2_data]
#N1_data = [float(i) for i in N1_data]
#N0_data = [float(i) for i in N0_data]

Q4_data = [float(i) for i in Q4_data]
Q3_data = [float(i) for i in Q3_data]
Q2_data = [float(i) for i in Q2_data]
Q1_data = [float(i) for i in Q1_data]
Q0_data = [float(i) for i in Q0_data]

mod_si_data = np.array(mod_si_data)
mod_b_data = np.array(mod_b_data)

#N3_data = np.array(N3_data)
N4_data = np.array(N4_data)
#N2_data = np.array(N2_data)
#N1_data = np.array(N1_data)
#N0_data = np.array(N0_data)

Q4_data = np.array(Q4_data)
Q3_data = np.array(Q3_data)
Q2_data = np.array(Q2_data)
Q1_data = np.array(Q1_data)
Q0_data = np.array(Q0_data)

mod_si_data_k2 = []
mod_b_data_k2 = []

#N3_data = []
N4_data_k2 = []
#N2_data = []
#N1_data = []
#N0_data = []

Q4_data_k2 = []
Q3_data_k2 = []
Q2_data_k2 = []
Q1_data_k2 = []
Q0_data_k2 = []



with open(f"{fil3}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        mod_si_data_k2.append(row[0])
        Q4_data_k2.append(row[1])
        Q3_data_k2.append(row[2])
        Q2_data_k2.append(row[3])
        Q1_data_k2.append(row[4])
        Q0_data_k2.append(row[5])

with open(f"{fil4}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        mod_b_data_k2.append(row[0])
#        N3_data.append(row[2])
        N4_data_k2.append(row[1])
#        N2_data.append(row[4])
#        N1_data.append(row[5])
#        N0_data.append(row[6])


mod_si_data_k2 = [float(i) for i in mod_si_data_k2]
mod_b_data_k2 = [float(i) for i in mod_b_data_k2]


#N3_data = [float(i) for i in N3_data]
N4_data_k2 = [float(i) for i in N4_data_k2]
#N2_data = [float(i) for i in N2_data]
#N1_data = [float(i) for i in N1_data]
#N0_data = [float(i) for i in N0_data]

Q4_data_k2 = [float(i) for i in Q4_data_k2]
Q3_data_k2 = [float(i) for i in Q3_data_k2]
Q2_data_k2 = [float(i) for i in Q2_data_k2]
Q1_data_k2 = [float(i) for i in Q1_data_k2]
Q0_data_k2 = [float(i) for i in Q0_data_k2]

mod_si_data_k2 = np.array(mod_si_data_k2)
mod_b_data_k2 = np.array(mod_b_data_k2)

#N3_data = np.array(N3_data)
N4_data_k2 = np.array(N4_data_k2)
#N2_data = np.array(N2_data)
#N1_data = np.array(N1_data)
#N0_data = np.array(N0_data)

Q4_data_k2 = np.array(Q4_data_k2)
Q3_data_k2 = np.array(Q3_data_k2)
Q2_data_k2 = np.array(Q2_data_k2)
Q1_data_k2 = np.array(Q1_data_k2)
Q0_data_k2 = np.array(Q0_data_k2)

mod_si_data_k05 = []
mod_b_data_k05 = []

#N3_data = []
N4_data_k05 = []
#N2_data = []
#N1_data = []
#N0_data = []

Q4_data_k05 = []
Q3_data_k05 = []
Q2_data_k05 = []
Q1_data_k05 = []
Q0_data_k05 = []



with open(f"{fil5}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        mod_si_data_k05.append(row[0])
        Q4_data_k05.append(row[1])
        Q3_data_k05.append(row[2])
        Q2_data_k05.append(row[3])
        Q1_data_k05.append(row[4])
        Q0_data_k05.append(row[5])

with open(f"{fil6}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        mod_b_data_k05.append(row[0])
#        N3_data.append(row[2])
        N4_data_k05.append(row[1])
#        N2_data.append(row[4])
#        N1_data.append(row[5])
#        N0_data.append(row[6])


mod_si_data_k05 = [float(i) for i in mod_si_data_k05]
mod_b_data_k05 = [float(i) for i in mod_b_data_k05]


#N3_data = [float(i) for i in N3_data]
N4_data_k05 = [float(i) for i in N4_data_k05]
#N2_data = [float(i) for i in N2_data]
#N1_data = [float(i) for i in N1_data]
#N0_data = [float(i) for i in N0_data]

Q4_data_k05 = [float(i) for i in Q4_data_k05]
Q3_data_k05 = [float(i) for i in Q3_data_k05]
Q2_data_k05 = [float(i) for i in Q2_data_k05]
Q1_data_k05 = [float(i) for i in Q1_data_k05]
Q0_data_k05 = [float(i) for i in Q0_data_k05]

mod_si_data_k05 = np.array(mod_si_data_k05)
mod_b_data_k05 = np.array(mod_b_data_k05)

#N3_data = np.array(N3_data)
N4_data_k05 = np.array(N4_data_k05)
#N2_data = np.array(N2_data)
#N1_data = np.array(N1_data)
#N0_data = np.array(N0_data)

Q4_data_k05 = np.array(Q4_data_k05)
Q3_data_k05 = np.array(Q3_data_k05)
Q2_data_k05 = np.array(Q2_data_k05)
Q1_data_k05 = np.array(Q1_data_k05)
Q0_data_k05 = np.array(Q0_data_k05)

mod_si_data_k6 = []
mod_b_data_k6 = []

#N3_data = []
N4_data_k6 = []
#N2_data = []
#N1_data = []
#N0_data = []

Q4_data_k6 = []
Q3_data_k6 = []
Q2_data_k6 = []
Q1_data_k6 = []
Q0_data_k6 = []



with open(f"{fil7}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        mod_si_data_k6.append(row[0])
        Q4_data_k6.append(row[1])
        Q3_data_k6.append(row[2])
        Q2_data_k6.append(row[3])
        Q1_data_k6.append(row[4])
        Q0_data_k6.append(row[5])

with open(f"{fil8}.csv", newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=';', quotechar='|')
    for row in spamreader:
        mod_b_data_k6.append(row[0])
#        N3_data.append(row[2])
        N4_data_k6.append(row[1])
#        N2_data.append(row[4])
#        N1_data.append(row[5])
#        N0_data.append(row[6])


mod_si_data_k6 = [float(i) for i in mod_si_data_k6]
mod_b_data_k6 = [float(i) for i in mod_b_data_k6]


#N3_data = [float(i) for i in N3_data]
N4_data_k6 = [float(i) for i in N4_data_k6]
#N2_data = [float(i) for i in N2_data]
#N1_data = [float(i) for i in N1_data]
#N0_data = [float(i) for i in N0_data]

Q4_data_k6 = [float(i) for i in Q4_data_k6]
Q3_data_k6 = [float(i) for i in Q3_data_k6]
Q2_data_k6 = [float(i) for i in Q2_data_k6]
Q1_data_k6 = [float(i) for i in Q1_data_k6]
Q0_data_k6 = [float(i) for i in Q0_data_k6]

mod_si_data_k6 = np.array(mod_si_data_k6)
mod_b_data_k6 = np.array(mod_b_data_k6)

#N3_data = np.array(N3_data)
N4_data_k6 = np.array(N4_data_k6)
#N2_data = np.array(N2_data)
#N1_data = np.array(N1_data)
#N0_data = np.array(N0_data)

Q4_data_k6 = np.array(Q4_data_k6)
Q3_data_k6 = np.array(Q3_data_k6)
Q2_data_k6 = np.array(Q2_data_k6)
Q1_data_k6 = np.array(Q1_data_k6)
Q0_data_k6 = np.array(Q0_data_k6)


#draw_nr = list(range(400))
#draw_ar = np.array(draw_nr)
#
#M2O = []
#
#for i in draw_ar:
#    next_mod = (((draw_ar[i]*0.5) / (100 + (draw_ar[i]*0.5)) * 100)* 2/3) + (((draw_ar[i]) / (100 + (draw_ar[i])) * 100) * 1/3)
#    M2O.append(next_mod)
#            #M1 = np.array((draw_ar/(100+draw_ar))*100)
#
##                w1 = [1, 0.01, 0.001]
#
#M2O.append(75)





def model(w1):
    
    
    H1 = [0, 8.4137,
          7.4205, 35.31424822, 
          0, 14.097, 
          22.898, 27.098]
    
    w_Tf = [math.exp(-(H1[0]/(0.008314462*w1[2]))), math.exp(-(H1[1]/(0.008314462*w1[2]))), 
            math.exp(-(H1[2]/(0.008314462*w1[2]))), math.exp(-(H1[3]/(0.008314462*w1[2]))), 
            math.exp(-(H1[4]/(0.008314462*w1[2]))), math.exp(-(H1[5]/(0.008314462*w1[2]))), 
            math.exp(-(H1[6]/(0.008314462*w1[2]))), math.exp(-(H1[7]/(0.008314462*w1[2])))] 
    
    w = np.array([w_Tf[0], w_Tf[1], w_Tf[2], w_Tf[3], w_Tf[4]*w1[0], w_Tf[5]*w1[0], w_Tf[6]*w1[0], w_Tf[7]*w1[0]])

    SSEv = []

    for i in list(range(len(k_si_data))):

        k_number = k_si_data[i]
        draw_range = (((1/(k_number+1))*300)+((k_number/(1+k_number))*400))
        draw_nr = list(range(int(draw_range)))
        draw_ar = np.array(draw_nr)

        M2O = []

        for i2 in draw_ar:
            next_mod = (((draw_ar[i2]) / ((k_number*0.5/(1+(k_number*0.5)))*200 + ((1 / (1 + (k_number*0.5)))*100) + draw_ar[i2]))*100)
            M2O.append(next_mod)

        B4_B2 = []
        for i3 in draw_ar:
            if M2O[i3] < w1[1]:
                B4_B2.append(1)
            else:
                B4_B2.append(0)
#((-((w1[7]*0.01)/((k_number*w1[7]*0.01)-k_number-1)))/((100)/(w[4]*k_number)))*100
        B4_B2 = np.array(B4_B2)

        B3 = [(1/((k_number*0.5)+1))*100, ]
        B4 = [0, ]
        B2 = [0, ]
        B1 = [0, ]
        B0 = [0, ]


        Q4 = [((k_number*0.5)/(1+(k_number*0.5)))*100, ]
        Q3 = [0, ]
        Q2 = [0, ]
        Q1 = [0, ]
        Q0 = [0, ]

        # The function that makes each iteration at a time


        for i4 in draw_ar:


            p_Q4 = Q4[-1]*w[4] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_Q3 = Q3[-1]*w[5] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_Q2 = Q2[-1]*w[6] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_Q1 = Q1[-1]*w[7] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

            p_B3 = B3[-1]*w[0] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_B4 = B4[-1]*w[1] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_B2 = B2[-1]*w[2] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_B1 = B1[-1]*w[3] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

            #Contribution to N4 from Si
            CQ4 = p_Q4 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CQ3 = p_Q3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CQ2 = p_Q2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CQ1 = p_Q1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

            #Contribution to N4 from B
            CB3 = p_B3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CB2 = p_B2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CB1 = p_B1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

            # Evolution of borate Qn units
            if B3[-1] - p_B3 - (p_B4 * CB3) < 0:
                next_B3 = 0
            else:
                next_B3 = B3[-1] - p_B3 - (p_B4 * CB3)

            if B4[-1] + (p_B3*B4_B2[i4]) - p_B4 < 0:
                next_B4 = 0
            else:
                next_B4 = B4[-1] + (p_B3*B4_B2[i4]) - p_B4

            if B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2) < 0:
                next_B2 = 0
            else:
                next_B2 = B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2)

            if B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
                next_B1 = 0
            else:
                next_B1 = B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

            if B0[-1] + p_B1 + (p_B4 * CB1 ) < 0:
                next_B0 = 0
            else:
                next_B0 = B0[-1] + p_B1 + (p_B4 * CB1)


            #Now for the evolution of silica Qn units

            if Q4[-1] - p_Q4 - (p_B4*CQ4) < 0:
                next_Q4 = 0
            else:
                next_Q4 = Q4[-1] - p_Q4 - (CQ4*p_B4)

            if Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4) < 0:
                next_Q3 = 0
            else:
                next_Q3 = Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4)

            if Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3) < 0:
                next_Q2 = 0
            else:
                next_Q2 = Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3)

            if Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2) < 0:
                next_Q1 = 0
            else:
                next_Q1 = Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2)

            if Q0[-1] + p_Q1 + (p_B4*CQ1) < 0:
                next_Q0 = 0
            else:
                next_Q0 = Q0[-1] + p_Q1 + (p_B4*CQ1)

            B3.append(next_B3)
            B4.append(next_B4)
            B2.append(next_B2)
            B1.append(next_B1)
            B0.append(next_B0)

            Q4.append(next_Q4)
            Q3.append(next_Q3)
            Q2.append(next_Q2)
            Q1.append(next_Q1)
            Q0.append(next_Q0)


        #Vi laver lister, som indeholer data fra modellen, som svarer til de punkter,
        #hvor vi har experimentiel data
        mod_si_m = []
#        B3_m = []
#        B4_m = []
#        B2_m = []
#        B1_m = []
#        B0_m = []

        Q4_m = []
        Q3_m = []
        Q2_m = []
        Q1_m = []
        Q0_m = []

        for i5 in mod_si_data:
            next_mod_m = min(M2O, key=lambda x:abs(x-i5))
            mod_si_m.append(next_mod_m)

#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B3_m = B3[ind]
#            B3_m.append(next_B3_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B4_m = B4[ind]
#            B4_m.append(next_B4_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B2_m = B2[ind]
#            B2_m.append(next_B2_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B1_m = B1[ind]
#            B1_m.append(next_B1_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B0_m = B0[ind]
#            B0_m.append(next_B0_m)




        ind = M2O.index(mod_si_m[i])
        next_Q4_m = Q4[ind]
        Q4_m.append(next_Q4_m)

        next_Q3_m = Q3[ind]
        Q3_m.append(next_Q3_m)

        next_Q2_m = Q2[ind]
        Q2_m.append(next_Q2_m)

        next_Q1_m = Q1[ind]
        Q1_m.append(next_Q1_m)

        next_Q0_m = Q0[ind]
        Q0_m.append(next_Q0_m)

#        B3_m = np.array(B3_m)
#        B4_m = np.array(B4_m)
#        B2_m = np.array(B2_m)
#        B1_m = np.array(B1_m)
#        B0_m = np.array(B0_m)

        Q4_m = np.array(Q4_m)
        Q3_m = np.array(Q3_m)
        Q2_m = np.array(Q2_m)
        Q1_m = np.array(Q1_m)
        Q0_m = np.array(Q0_m)


        SSEp = sum(((Q4_data[i] - Q4_m)**2) + ((Q3_data[i] - Q3_m)**2) + ((Q2_data[i] - Q2_m)**2) + ((Q1_data[i] - Q1_m)**2) + ((Q0_data[i] - Q0_m)**2))
        SSEv.append(SSEp)


    for i in list(range(len(k_b_data))):

        k_number = k_b_data[i]
        draw_range = (((1/(k_number+1))*300)+((k_number/(1+k_number))*400))
        draw_nr = list(range(int(draw_range)))
        draw_ar = np.array(draw_nr)

        M2O = []

        for i2 in draw_ar:
            next_mod = (((draw_ar[i2]) / ((k_number*0.5/(1+(k_number*0.5)))*200 + ((1 / (1 + (k_number*0.5)))*100) + draw_ar[i2]))*100)
            M2O.append(next_mod)

        B4_B2 = []
        for i3 in draw_ar:
            if M2O[i3] < w1[1]:
                B4_B2.append(1)
            else:
                B4_B2.append(0)

        B4_B2 = np.array(B4_B2)

        B3 = [(1/((k_number*0.5)+1))*100, ]
        B4 = [0, ]
        B2 = [0, ]
        B1 = [0, ]
        B0 = [0, ]


        Q4 = [((k_number*0.5)/(1+(k_number*0.5)))*100, ]
        Q3 = [0, ]
        Q2 = [0, ]
        Q1 = [0, ]
        Q0 = [0, ]

        # The function that makes each iteration at a time


        for i4 in draw_ar:


            p_Q4 = Q4[-1]*w[4] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_Q3 = Q3[-1]*w[5] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_Q2 = Q2[-1]*w[6] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_Q1 = Q1[-1]*w[7] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

            p_B3 = B3[-1]*w[0] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_B4 = B4[-1]*w[1] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_B2 = B2[-1]*w[2] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
            p_B1 = B1[-1]*w[3] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

            #Contribution to N4 from Si
            CQ4 = p_Q4 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CQ3 = p_Q3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CQ2 = p_Q2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CQ1 = p_Q1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

            #Contribution to N4 from B
            CB3 = p_B3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CB2 = p_B2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
            CB1 = p_B1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

            # Evolution of borate Qn units
            if B3[-1] - p_B3 - (p_B4 * CB3) < 0:
                next_B3 = 0
            else:
                next_B3 = B3[-1] - p_B3 - (p_B4 * CB3)

            if B4[-1] + (p_B3*B4_B2[i4]) - p_B4 < 0:
                next_B4 = 0
            else:
                next_B4 = B4[-1] + (p_B3*B4_B2[i4]) - p_B4

            if B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2) < 0:
                next_B2 = 0
            else:
                next_B2 = B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2)

            if B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
                next_B1 = 0
            else:
                next_B1 = B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

            if B0[-1] + p_B1 + (p_B4 * CB1 ) < 0:
                next_B0 = 0
            else:
                next_B0 = B0[-1] + p_B1 + (p_B4 * CB1)


            #Now for the evolution of silica Qn units

            if Q4[-1] - p_Q4 - (p_B4*CQ4) < 0:
                next_Q4 = 0
            else:
                next_Q4 = Q4[-1] - p_Q4 - (CQ4*p_B4)

            if Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4) < 0:
                next_Q3 = 0
            else:
                next_Q3 = Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4)

            if Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3) < 0:
                next_Q2 = 0
            else:
                next_Q2 = Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3)

            if Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2) < 0:
                next_Q1 = 0
            else:
                next_Q1 = Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2)

            if Q0[-1] + p_Q1 + (p_B4*CQ1) < 0:
                next_Q0 = 0
            else:
                next_Q0 = Q0[-1] + p_Q1 + (p_B4*CQ1)

            B3.append(next_B3)
            B4.append(next_B4)
            B2.append(next_B2)
            B1.append(next_B1)
            B0.append(next_B0)

            Q4.append(next_Q4)
            Q3.append(next_Q3)
            Q2.append(next_Q2)
            Q1.append(next_Q1)
            Q0.append(next_Q0)


        #Vi laver lister, som indeholer data fra modellen, som svarer til de punkter,
        #hvor vi har experimentiel data
        mod_b_m = []
#        B3_m = []
        B4_m = []
#        B2_m = []
#        B1_m = []
#        B0_m = []

#        Q4_m = []
#        Q3_m = []
#        Q2_m = []
#        Q1_m = []
#        Q0_m = []

        for i5 in mod_b_data:
            next_mod_m = min(M2O, key=lambda x:abs(x-i5))
            mod_b_m.append(next_mod_m)

#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B3_m = B3[ind]
#            B3_m.append(next_B3_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B4_m = B4[ind]
#            B4_m.append(next_B4_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B2_m = B2[ind]
#            B2_m.append(next_B2_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B1_m = B1[ind]
#            B1_m.append(next_B1_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B0_m = B0[ind]
#            B0_m.append(next_B0_m)




        ind = M2O.index(mod_b_m[i])
#        next_B3_m = B3[ind]
#        B3_m.append(next_B3_m)

        next_B4_m = B4[ind]
        B4_m.append(next_B4_m)

#        next_B2_m = B2[ind]
#        B2_m.append(next_B2_m)
#
#        next_B1_m = B1[ind]
#        B1_m.append(next_B1_m)
#
#        next_B0_m = B0[ind]
#        B0_m.append(next_B0_m)

#        B3_m = np.array(B3_m)
        B4_m = np.array(B4_m)
#        B2_m = np.array(B2_m)
#        B1_m = np.array(B1_m)
#        B0_m = np.array(B0_m)

#        Q4_m = np.array(Q4_m)
#        Q3_m = np.array(Q3_m)
#        Q2_m = np.array(Q2_m)
#        Q1_m = np.array(Q1_m)
#        Q0_m = np.array(Q0_m)


        SSEs = sum(((N4_data[i] - B4_m)**2))*5
        SSEv.append(SSEs)



    SSE = sum(SSEv)

    return SSE

#res = scipy.optimize.minimize(model, w0, method='nelder-mead', options={'xtol': 1e-8, 'disp': True})

#Hvordan inkluderer jeg, at w[1] er K afhængig? uden at kende afhængigheden
w0 = [0.163307865824085, 21.57506806662955, 399.492879273407]
w0 = np.array(w0)

#class MyBounds(object):
#    def __init__(self, xmax=[100000., 100000., 100000., 100000., 100000., 100000., 100000., 100.], xmin=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ):
#        self.xmax = np.array(xmax)
#        self.xmin = np.array(xmin)
#    def __call__(self, **kwargs):
#        x = kwargs["x_new"]
#        tmax = bool(np.all(x <= self.xmax))
#        tmin = bool(np.all(x >= self.xmin))
#        return tmax and tmin


minimizer_kwargs = {"method": "BFGS"}

#mybounds = MyBounds()

SSElist =[]
it = 1

def print_fun(x, f, accepted):
    global it
    global SSElist
    print("SSE = {} for {}".format(f, it))
    it += 1
    SSElist.append(f)

res = scipy.optimize.basinhopping(model, w0, niter=1, T=1.0, stepsize=0.2, minimizer_kwargs=minimizer_kwargs, take_step=None, accept_test=None, callback=print_fun, interval=50, disp=True, niter_success=None, seed=None)

print(f"Fitte parametre NaBSi: wsi/b = {res.x[0]}, B4/B2 = {res.x[1]}, Tf = {res.x[2]}")

H1 = [0, 8.4137,
      7.4205, 35.31424822, 
      0, 14.097, 
      22.898, 27.098]

w_Tf = [math.exp(-(H1[0]/(0.008314462*res.x[2]))), math.exp(-(H1[1]/(0.008314462*res.x[2]))), 
        math.exp(-(H1[2]/(0.008314462*res.x[2]))), math.exp(-(H1[3]/(0.008314462*res.x[2]))), 
        math.exp(-(H1[4]/(0.008314462*res.x[2]))), math.exp(-(H1[5]/(0.008314462*res.x[2]))), 
        math.exp(-(H1[6]/(0.008314462*res.x[2]))), math.exp(-(H1[7]/(0.008314462*res.x[2])))] 

w = np.array([w_Tf[0], w_Tf[1], w_Tf[2], w_Tf[3], w_Tf[4]*res.x[0], w_Tf[5]*res.x[0], w_Tf[6]*res.x[0], w_Tf[7]*res.x[0]])

SSEv = []

Q4_mp = []
Q3_mp = []
Q2_mp = []
Q1_mp = []
Q0_mp = []

N4_mp = []

for i in list(range(len(k_si_data))):

    k_number = k_si_data[i]
    draw_range = (((1/(k_number+1))*300)+((k_number/(1+k_number))*400))
    draw_nr = list(range(int(draw_range)))
    draw_ar = np.array(draw_nr)

    M2O = []

    for i2 in draw_ar:
        next_mod = (((draw_ar[i2]) / ((k_number*0.5/(1+(k_number*0.5)))*200 + ((1 / (1 + (k_number*0.5)))*100) + draw_ar[i2]))*100)
        M2O.append(next_mod)

    B4_B2 = []
    for i3 in draw_ar:
        if M2O[i3] < res.x[1]:
            B4_B2.append(1)
        else:
            B4_B2.append(0)

    B4_B2 = np.array(B4_B2)

    B3 = [(1/((k_number*0.5)+1))*100, ]
    B4 = [0, ]
    B2 = [0, ]
    B1 = [0, ]
    B0 = [0, ]


    Q4 = [((k_number*0.5)/(1+(k_number*0.5)))*100, ]
    Q3 = [0, ]
    Q2 = [0, ]
    Q1 = [0, ]
    Q0 = [0, ]

    # The function that makes each iteration at a time


    for i4 in draw_ar:


        p_Q4 = Q4[-1]*w[4] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_Q3 = Q3[-1]*w[5] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_Q2 = Q2[-1]*w[6] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_Q1 = Q1[-1]*w[7] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

        p_B3 = B3[-1]*w[0] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_B4 = B4[-1]*w[1] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_B2 = B2[-1]*w[2] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_B1 = B1[-1]*w[3] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

        #Contribution to N4 from Si
        CQ4 = p_Q4 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CQ3 = p_Q3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CQ2 = p_Q2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CQ1 = p_Q1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

        #Contribution to N4 from B
        CB3 = p_B3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CB2 = p_B2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CB1 = p_B1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

        # Evolution of borate Qn units
        if B3[-1] - p_B3 - (p_B4 * CB3) < 0:
            next_B3 = 0
        else:
            next_B3 = B3[-1] - p_B3 - (p_B4 * CB3)

        if B4[-1] + (p_B3*B4_B2[i4]) - p_B4 < 0:
            next_B4 = 0
        else:
            next_B4 = B4[-1] + (p_B3*B4_B2[i4]) - p_B4

        if B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2) < 0:
            next_B2 = 0
        else:
            next_B2 = B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2)

        if B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
            next_B1 = 0
        else:
            next_B1 = B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

        if B0[-1] + p_B1 + (p_B4 * CB1 ) < 0:
            next_B0 = 0
        else:
            next_B0 = B0[-1] + p_B1 + (p_B4 * CB1)


        #Now for the evolution of silica Qn units

        if Q4[-1] - p_Q4 - (p_B4*CQ4) < 0:
            next_Q4 = 0
        else:
            next_Q4 = Q4[-1] - p_Q4 - (CQ4*p_B4)

        if Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4) < 0:
            next_Q3 = 0
        else:
            next_Q3 = Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4)

        if Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3) < 0:
            next_Q2 = 0
        else:
            next_Q2 = Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3)

        if Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2) < 0:
            next_Q1 = 0
        else:
            next_Q1 = Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2)

        if Q0[-1] + p_Q1 + (p_B4*CQ1) < 0:
            next_Q0 = 0
        else:
            next_Q0 = Q0[-1] + p_Q1 + (p_B4*CQ1)

        B3.append(next_B3)
        B4.append(next_B4)
        B2.append(next_B2)
        B1.append(next_B1)
        B0.append(next_B0)

        Q4.append(next_Q4)
        Q3.append(next_Q3)
        Q2.append(next_Q2)
        Q1.append(next_Q1)
        Q0.append(next_Q0)


    #Vi laver lister, som indeholer data fra modellen, som svarer til de punkter,
    #hvor vi har experimentiel data
    mod_si_m = []
#        B3_m = []
#        B4_m = []
#        B2_m = []
#        B1_m = []
#        B0_m = []



    for i5 in mod_si_data:
        next_mod_m = min(M2O, key=lambda x:abs(x-i5))
        mod_si_m.append(next_mod_m)

#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B3_m = B3[ind]
#            B3_m.append(next_B3_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B4_m = B4[ind]
#            B4_m.append(next_B4_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B2_m = B2[ind]
#            B2_m.append(next_B2_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B1_m = B1[ind]
#            B1_m.append(next_B1_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B0_m = B0[ind]
#            B0_m.append(next_B0_m)




    ind = M2O.index(mod_si_m[i])
    next_Q4_m = Q4[ind]
    Q4_mp.append(next_Q4_m)

    next_Q3_m = Q3[ind]
    Q3_mp.append(next_Q3_m)

    next_Q2_m = Q2[ind]
    Q2_mp.append(next_Q2_m)

    next_Q1_m = Q1[ind]
    Q1_mp.append(next_Q1_m)

    next_Q0_m = Q0[ind]
    Q0_mp.append(next_Q0_m)

#        B3_m = np.array(B3_m)
#        B4_m = np.array(B4_m)
#        B2_m = np.array(B2_m)
#        B1_m = np.array(B1_m)
#        B0_m = np.array(B0_m)





for i in list(range(len(k_b_data))):

    k_number = k_b_data[i]
    draw_range = (((1/(k_number+1))*300)+((k_number/(1+k_number))*400))
    draw_nr = list(range(int(draw_range)))
    draw_ar = np.array(draw_nr)

    M2O = []

    for i2 in draw_ar:
        next_mod = (((draw_ar[i2]) / ((k_number*0.5/(1+(k_number*0.5)))*200 + ((1 / (1 + (k_number*0.5)))*100) + draw_ar[i2]))*100)
        M2O.append(next_mod)

    B4_B2 = []
    for i3 in draw_ar:
        if M2O[i3] < res.x[1]:
            B4_B2.append(1)
        else:
            B4_B2.append(0)

    B4_B2 = np.array(B4_B2)

    B3 = [(1/((k_number*0.5)+1))*100, ]
    B4 = [0, ]
    B2 = [0, ]
    B1 = [0, ]
    B0 = [0, ]


    Q4 = [((k_number*0.5)/(1+(k_number*0.5)))*100, ]
    Q3 = [0, ]
    Q2 = [0, ]
    Q1 = [0, ]
    Q0 = [0, ]

    # The function that makes each iteration at a time


    for i4 in draw_ar:


        p_Q4 = Q4[-1]*w[4] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_Q3 = Q3[-1]*w[5] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_Q2 = Q2[-1]*w[6] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_Q1 = Q1[-1]*w[7] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

        p_B3 = B3[-1]*w[0] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_B4 = B4[-1]*w[1] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_B2 = B2[-1]*w[2] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])
        p_B1 = B1[-1]*w[3] / (Q4[-1]*w[4] + Q3[-1]*w[5] + Q2[-1]*w[6] + Q1[-1]*w[7] + B3[-1]*w[0] + B4[-1]*w[1] + B2[-1]*w[2] + B1[-1]*w[3])

        #Contribution to N4 from Si
        CQ4 = p_Q4 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CQ3 = p_Q3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CQ2 = p_Q2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CQ1 = p_Q1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

        #Contribution to N4 from B
        CB3 = p_B3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CB2 = p_B2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
        CB1 = p_B1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

        # Evolution of borate Qn units
        if B3[-1] - p_B3 - (p_B4 * CB3) < 0:
            next_B3 = 0
        else:
            next_B3 = B3[-1] - p_B3 - (p_B4 * CB3)

        if B4[-1] + (p_B3*B4_B2[i4]) - p_B4 < 0:
            next_B4 = 0
        else:
            next_B4 = B4[-1] + (p_B3*B4_B2[i4]) - p_B4

        if B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2) < 0:
            next_B2 = 0
        else:
            next_B2 = B2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2)

        if B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
            next_B1 = 0
        else:
            next_B1 = B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

        if B0[-1] + p_B1 + (p_B4 * CB1 ) < 0:
            next_B0 = 0
        else:
            next_B0 = B0[-1] + p_B1 + (p_B4 * CB1)


        #Now for the evolution of silica Qn units

        if Q4[-1] - p_Q4 - (p_B4*CQ4) < 0:
            next_Q4 = 0
        else:
            next_Q4 = Q4[-1] - p_Q4 - (CQ4*p_B4)

        if Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4) < 0:
            next_Q3 = 0
        else:
            next_Q3 = Q3[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4)

        if Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3) < 0:
            next_Q2 = 0
        else:
            next_Q2 = Q2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3)

        if Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2) < 0:
            next_Q1 = 0
        else:
            next_Q1 = Q1[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2)

        if Q0[-1] + p_Q1 + (p_B4*CQ1) < 0:
            next_Q0 = 0
        else:
            next_Q0 = Q0[-1] + p_Q1 + (p_B4*CQ1)

        B3.append(next_B3)
        B4.append(next_B4)
        B2.append(next_B2)
        B1.append(next_B1)
        B0.append(next_B0)

        Q4.append(next_Q4)
        Q3.append(next_Q3)
        Q2.append(next_Q2)
        Q1.append(next_Q1)
        Q0.append(next_Q0)


    #Vi laver lister, som indeholer data fra modellen, som svarer til de punkter,
    #hvor vi har experimentiel data
    mod_b_m = []
#        B3_m = []
#    B4_m = []
#        B2_m = []
#        B1_m = []
#        B0_m = []

#        Q4_m = []
#        Q3_m = []
#        Q2_m = []
#        Q1_m = []
#        Q0_m = []

    for i5 in mod_b_data:
        next_mod_m = min(M2O, key=lambda x:abs(x-i5))
        mod_b_m.append(next_mod_m)

#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B3_m = B3[ind]
#            B3_m.append(next_B3_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B4_m = B4[ind]
#            B4_m.append(next_B4_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B2_m = B2[ind]
#            B2_m.append(next_B2_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B1_m = B1[ind]
#            B1_m.append(next_B1_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B0_m = B0[ind]
#            B0_m.append(next_B0_m)




    ind = M2O.index(mod_b_m[i])
#        next_B3_m = B3[ind]
#        B3_m.append(next_B3_m)

    next_B4_m = B4[ind]
    N4_mp.append(next_B4_m)

#        next_B2_m = B2[ind]
#        B2_m.append(next_B2_m)
#
#        next_B1_m = B1[ind]
#        B1_m.append(next_B1_m)
#
#        next_B0_m = B0[ind]
#        B0_m.append(next_B0_m)

#        B3_m = np.array(B3_m)

#        B2_m = np.array(B2_m)
#        B1_m = np.array(B1_m)
#        B0_m = np.array(B0_m)

#        Q4_m = np.array(Q4_m)
#        Q3_m = np.array(Q3_m)
#        Q2_m = np.array(Q2_m)
#        Q1_m = np.array(Q1_m)
#        Q0_m = np.array(Q0_m)



Q4_mp = np.array(Q4_mp)
Q3_mp = np.array(Q3_mp)
Q2_mp = np.array(Q2_mp)
Q1_mp = np.array(Q1_mp)
Q0_mp = np.array(Q0_mp)

N4_mp = np.array(N4_mp)

SSEp = sum(((Q4_data - Q4_mp)**2) + ((Q3_data - Q3_mp)**2) + ((Q2_data - Q2_mp)**2) + ((Q1_data - Q1_mp)**2) + ((Q0_data - Q0_mp)**2))
SSEs = sum(((N4_data - N4_mp)**2))
SSEv.append(SSEs)
SSEv.append(SSEp)

SSE = sum(SSEv)
t = list(range(100))
plt.plot(Q4_data, Q4_mp, 'rd', Q3_data, Q3_mp, 'bd', Q2_data, Q2_mp, 'gd', Q1_data, Q1_mp, 'yd', N4_data, N4_mp, 'kd', 
#         Q0_data, Q0_mp, 'cd'
         )
plt.plot(t, t, 'k--')
plt.axis([0, 65, 0, 65])
plt.xlabel("Experimental Data")
plt.ylabel("Model Data")
plt.title('Quality of fit')
plt.savefig(os.path.join(f'{fil}_data',f'{fil} Borosilicate_1w_Tf.png'))
plt.show()
#
print(f"Parameters: ")
print(f"SSE: {SSE}")


        
m_data_si_plt_kval = np.column_stack([Q4_data, Q4_mp, Q3_data, Q3_mp, Q2_data, Q2_mp, Q1_data, Q1_mp])
m_data_b_plt_kval = np.column_stack([N4_data, N4_mp])
np.savetxt(os.path.join(f'{fil}_data', f"Model_Qn_v_data_Qn_{fil}.csv"), m_data_si_plt_kval)
np.savetxt(os.path.join(f'{fil}_data', f"Model_Bn_v_data_Qn_{fil}.csv"), m_data_b_plt_kval)
#
#plt.plot(M2O, Q4, 'r-', M2O, Q3, 'b-', M2O, Q2, 'g-', M2O,Q1, 'y-', M2O, Q0, 'k-')
#plt.plot(mod_data, Q4_data, 'rd', mod_data, Q3_data, 'bd', mod_data, Q2_data, 'gd', mod_data, Q1_data, 'yd', mod_data, Q0_data, 'kd')
#plt.axis([0, 75, 0, 40])
#plt.xlabel("Modifier mol %")
#plt.ylabel(f"Qn species concentration")
#plt.title('Qn distribution')
#plt.savefig(f'{fil} Silicate.png')
#plt.show()
#
##M2O = np.array(M2O)
#
#np.savetxt(f"params {fil}_Tf", res.x)
#
##m_data = np.column_stack([M2O, B3, T3, T4, D, L4, M3, P3, O3])
##np.savetxt(f"Model data {fil}.csv", m_data)
#
##m_Qdata = np.column_stack([M2O, Q3, Q4, Q2, Q1, Q0])
##np.savetxt(f"Model Qdata {fil}.csv", m_data)
#

#k_si_set = list(set(k_si_data))
#k_b_set = list(set(k_b_data))
#
#
#
#k_si_ind = [i for i, n in enumerate(k_si_data) if n == 0.5]

#Q4_mp_k2 = []
#Q3_mp_k2 = []
#Q2_mp_k2 = []
#Q1_mp_k2 = []
#Q0_mp_k2 = []
#
#N4_mp_k2 = []

k_number = 2
draw_range = (((1/(k_number+1))*300)+((k_number/(1+k_number))*400))
draw_nr = list(range(int(draw_range)))
draw_ar = np.array(draw_nr)

M2O_k2 = []

for i2 in draw_ar:
    next_mod = (((draw_ar[i2]) / ((k_number*0.5/(1+(k_number*0.5)))*200 + ((1 / (1 + (k_number*0.5)))*100) + draw_ar[i2]))*100)
    M2O_k2.append(next_mod)

M2O_k2.append(75)

B4_B2 = []
for i3 in draw_ar:
    if M2O_k2[i3] < res.x[1]:
        B4_B2.append(1)
    else:
        B4_B2.append(0)


B4_B2 = np.array(B4_B2)

B3_k2 = [(1/((k_number*0.5)+1))*100, ]
B4_k2 = [0, ]
B2_k2 = [0, ]
B1_k2 = [0, ]
B0_k2 = [0, ]


Q4_k2 = [((k_number*0.5)/(1+(k_number*0.5)))*100, ]
Q3_k2 = [0, ]
Q2_k2 = [0, ]
Q1_k2 = [0, ]
Q0_k2 = [0, ]

# The function that makes each iteration at a time


for i4 in draw_ar:


    p_Q4 = Q4_k2[-1]*w[4] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])
    p_Q3 = Q3_k2[-1]*w[5] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])
    p_Q2 = Q2_k2[-1]*w[6] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])
    p_Q1 = Q1_k2[-1]*w[7] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])
    
    p_B3 = B3_k2[-1]*w[0] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])
    p_B4 = B4_k2[-1]*w[1] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])
    p_B2 = B2_k2[-1]*w[2] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])
    p_B1 = B1_k2[-1]*w[3] / (Q4_k2[-1]*w[4] + Q3_k2[-1]*w[5] + Q2_k2[-1]*w[6] + Q1_k2[-1]*w[7] + B3_k2[-1]*w[0] + B4_k2[-1]*w[1] + B2_k2[-1]*w[2] + B1_k2[-1]*w[3])

    #Contribution to N4 from Si
    CQ4 = p_Q4 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ3 = p_Q3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ2 = p_Q2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ1 = p_Q1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

    #Contribution to N4 from B
    CB3 = p_B3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CB2 = p_B2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CB1 = p_B1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

    # Evolution of borate Qn units
    if B3_k2[-1] - p_B3 - (p_B4 * CB3) < 0:
        next_B3 = 0
    else:
        next_B3 = B3_k2[-1] - p_B3 - (p_B4 * CB3)

    if B4_k2[-1] + (p_B3*B4_B2[i4]) - p_B4 < 0:
        next_B4 = 0
    else:
        next_B4 = B4_k2[-1] + (p_B3*B4_B2[i4]) - p_B4

    if B2_k2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2) < 0:
        next_B2 = 0
    else:
        next_B2 = B2_k2[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2)

    if B1_k2[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
        next_B1 = 0
    else:
        next_B1 = B1_k2[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

    if B0_k2[-1] + p_B1 + (p_B4 * CB1 ) < 0:
        next_B0 = 0
    else:
        next_B0 = B0_k2[-1] + p_B1 + (p_B4 * CB1)


    #Now for the evolution of silica Qn units

    if Q4_k2[-1] - p_Q4 - (p_B4*CQ4) < 0:
        next_Q4 = 0
    else:
        next_Q4 = Q4_k2[-1] - p_Q4 - (CQ4*p_B4)

    if Q3_k2[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4) < 0:
        next_Q3 = 0
    else:
        next_Q3 = Q3_k2[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4)

    if Q2_k2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3) < 0:
        next_Q2 = 0
    else:
        next_Q2 = Q2_k2[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3)

    if Q1_k2[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2) < 0:
        next_Q1 = 0
    else:
        next_Q1 = Q1_k2[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2)

    if Q0_k2[-1] + p_Q1 + (p_B4*CQ1) < 0:
        next_Q0 = 0
    else:
        next_Q0 = Q0_k2[-1] + p_Q1 + (p_B4*CQ1)

    B3_k2.append(next_B3)
    B4_k2.append(next_B4)
    B2_k2.append(next_B2)
    B1_k2.append(next_B1)
    B0_k2.append(next_B0)

    Q4_k2.append(next_Q4)
    Q3_k2.append(next_Q3)
    Q2_k2.append(next_Q2)
    Q1_k2.append(next_Q1)
    Q0_k2.append(next_Q0)


#Vi laver lister, som indeholer data fra modellen, som svarer til de punkter,
#hvor vi har experimentiel data
#mod_si_m_k2 = []
##        B3_m = []
##        B4_m = []
##        B2_m = []
##        B1_m = []
##        B0_m = []
#
#
#
#for i5 in mod_si_data_k2:
#    next_mod_m = min(M2O_k2, key=lambda x:abs(x-i5))
#    mod_si_m_k2.append(next_mod_m)

#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B3_m = B3[ind]
#            B3_m.append(next_B3_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B4_m = B4[ind]
#            B4_m.append(next_B4_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B2_m = B2[ind]
#            B2_m.append(next_B2_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B1_m = B1[ind]
#            B1_m.append(next_B1_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B0_m = B0[ind]
#            B0_m.append(next_B0_m)

#==============================================================================
#KOMMET HER TIL!!!
#==============================================================================

#for i in range(mod_si_m_k2):
#    ind = M2O_k2.index(mod_si_m_k2[i])
#    next_Q4_m = Q4_k2[ind]
#    Q4_mp_k2.append(next_Q4_m)
#
#    next_Q3_m = Q3_k2[ind]
#    Q3_mp_k2.append(next_Q3_m)
#
#    next_Q2_m = Q2_k2[ind]
#    Q2_mp_k2.append(next_Q2_m)
#
#    next_Q1_m = Q1_k2[ind]
#    Q1_mp_k2.append(next_Q1_m)
#
#    next_Q0_m = Q0_k2[ind]
#    Q0_mp_k2.append(next_Q0_m)

#        B3_m = np.array(B3_m)
#        B4_m = np.array(B4_m)
#        B2_m = np.array(B2_m)
#        B1_m = np.array(B1_m)
#        B0_m = np.array(B0_m)

#mod_b_m_k2 = []
#
#
#for i5 in mod_b_data_k2:
#    next_mod_m = min(M2O_k2, key=lambda x:abs(x-i5))
#    mod_b_m_k2.append(next_mod_m)

#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B3_m = B3[ind]
#            B3_m.append(next_B3_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B4_m = B4[ind]
#            B4_m.append(next_B4_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B2_m = B2[ind]
#            B2_m.append(next_B2_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B1_m = B1[ind]
#            B1_m.append(next_B1_m)
#
#        for i in mod_m:
#            ind = M2O.index(i)
#            next_B0_m = B0[ind]
#            B0_m.append(next_B0_m)




#for i in range(mod_b_m_k2):
#    ind = M2O_k2.index(mod_b_m_k2[i])
#    next_B4_m = Q4_k2[ind]
#    N4_mp_k2.append(next_B4_m)
#        next_B2_m = B2[ind]
#        B2_m.append(next_B2_m)
#
#        next_B1_m = B1[ind]
#        B1_m.append(next_B1_m)
#
#        next_B0_m = B0[ind]
#        B0_m.append(next_B0_m)

#        B3_m = np.array(B3_m)

#        B2_m = np.array(B2_m)
#        B1_m = np.array(B1_m)
#        B0_m = np.array(B0_m)

#        Q4_m = np.array(Q4_m)
#        Q3_m = np.array(Q3_m)
#        Q2_m = np.array(Q2_m)
#        Q1_m = np.array(Q1_m)
#        Q0_m = np.array(Q0_m)
    

plt.plot(M2O_k2, Q4_k2, 'r-', M2O_k2, Q3_k2, 'b-', M2O_k2, Q2_k2, 'g-', M2O_k2, Q1_k2, 'y-', M2O_k2, Q0_k2, 'k-')
plt.plot(mod_si_data_k2, Q4_data_k2, 'rd', mod_si_data_k2, Q3_data_k2, 'bd', mod_si_data_k2, Q2_data_k2, 'gd', mod_si_data_k2, Q1_data_k2, 'yd', mod_si_data_k2, Q0_data_k2, 'kd')
plt.axis([0, 65, 0, 65])
plt.xlabel("Modifier mol %")
plt.ylabel(f"Qn species concentration")
plt.title('Qn distribution K = 2')
plt.savefig(os.path.join(f'{fil}_data',f'{fil} Silicate_k2_0w.png'))
plt.show()

plt.plot(
#        M2O_k2, B3_k2, 'r-', 
        M2O_k2, B4_k2, 'b-'
#        , M2O_k2, B2_k2, 'g-', M2O_k2, B1_k2, 'y-', M2O_k2, B0_k2, 'k-'
        )
plt.plot(mod_b_data_k2, N4_data_k2, 'bd')
plt.axis([0, 65, 0, 65])
plt.xlabel("Modifier mol %")
plt.ylabel(f"Bn species concentration")
plt.title('Bn distribution K = 2')
plt.savefig(os.path.join(f'{fil}_data', f'{fil} Borate_k2_0w.png'))
plt.show()


        
m_Qmod_k2 = np.column_stack([M2O_k2, Q4_k2, Q3_k2, Q2_k2, Q1_k2, Q0_k2])
np.savetxt(os.path.join(f'{fil}_data', f"k2_Qn__model_{fil}.csv"), m_Qmod_k2)

m_Qdata_k2 = np.column_stack([mod_si_data_k2, Q4_data_k2, Q3_data_k2, Q2_data_k2, Q1_data_k2, Q0_data_k2])
np.savetxt(os.path.join(f'{fil}_data', f"k2_Qn__data_{fil}.csv"), m_Qdata_k2)

m_Bmod_k2 = np.column_stack([M2O_k2, B4_k2])
np.savetxt(os.path.join(f'{fil}_data', f"k2_Bn_model_{fil}.csv"), m_Bmod_k2)

m_Bdata_k2 = np.column_stack([mod_b_data_k2, N4_data_k2])
np.savetxt(os.path.join(f'{fil}_data', f"k2_Bn_data_{fil}.csv"), m_Bdata_k2)


##m_Qdata = np.column_stack([M2O, Q3, Q4, Q2, Q1, Q0])
##np.savetxt(f"Model Qdata {fil}.csv", m_data)



Q4_mp_k05 = []
Q3_mp_k05 = []
Q2_mp_k05 = []
Q1_mp_k05 = []
Q0_mp_k05 = []

N4_mp_k05 = []

k_number = 4
draw_range = (((1/(k_number+1))*300)+((k_number/(1+k_number))*400))
draw_nr = list(range(int(draw_range)))
draw_ar = np.array(draw_nr)

M2O_k05 = []

for i2 in draw_ar:
    next_mod = (((draw_ar[i2]) / ((k_number*0.5/(1+(k_number*0.5)))*200 + ((1 / (1 + (k_number*0.5)))*100) + draw_ar[i2]))*100)
    M2O_k05.append(next_mod)

M2O_k05.append(75)

B4_B2 = []
for i3 in draw_ar:
    if M2O_k05[i3] < res.x[1]:
        B4_B2.append(1)
    else:
        B4_B2.append(0)

B4_B2 = np.array(B4_B2)

B3_k05 = [(1/((k_number*0.5)+1))*100, ]
B4_k05 = [0, ]
B2_k05 = [0, ]
B1_k05 = [0, ]
B0_k05 = [0, ]


Q4_k05 = [((k_number*0.5)/(1+(k_number*0.5)))*100, ]
Q3_k05 = [0, ]
Q2_k05 = [0, ]
Q1_k05 = [0, ]
Q0_k05 = [0, ]

# The function that makes each iteration at a time


for i4 in draw_ar:


    p_Q4 = Q4_k05[-1]*w[4] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])
    p_Q3 = Q3_k05[-1]*w[5] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])
    p_Q2 = Q2_k05[-1]*w[6] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])
    p_Q1 = Q1_k05[-1]*w[7] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])
    
    p_B3 = B3_k05[-1]*w[0] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])
    p_B4 = B4_k05[-1]*w[1] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])
    p_B2 = B2_k05[-1]*w[2] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])
    p_B1 = B1_k05[-1]*w[3] / (Q4_k05[-1]*w[4] + Q3_k05[-1]*w[5] + Q2_k05[-1]*w[6] + Q1_k05[-1]*w[7] + B3_k05[-1]*w[0] + B4_k05[-1]*w[1] + B2_k05[-1]*w[2] + B1_k05[-1]*w[3])

    #Contribution to N4 from Si
    CQ4 = p_Q4 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ3 = p_Q3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ2 = p_Q2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ1 = p_Q1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

    #Contribution to N4 from B
    CB3 = p_B3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CB2 = p_B2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CB1 = p_B1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

    # Evolution of borate Qn units
    if B3_k05[-1] - p_B3 - (p_B4 * CB3) < 0:
        next_B3 = 0
    else:
        next_B3 = B3_k05[-1] - p_B3 - (p_B4 * CB3)

    if B4_k05[-1] + (p_B3*B4_B2[i4]) - p_B4 < 0:
        next_B4 = 0
    else:
        next_B4 = B4_k05[-1] + (p_B3*B4_B2[i4]) - p_B4

    if B2_k05[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2) < 0:
        next_B2 = 0
    else:
        next_B2 = B2_k05[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2)

    if B1_k05[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
        next_B1 = 0
    else:
        next_B1 = B1_k05[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

    if B0_k05[-1] + p_B1 + (p_B4 * CB1 ) < 0:
        next_B0 = 0
    else:
        next_B0 = B0_k05[-1] + p_B1 + (p_B4 * CB1)


    #Now for the evolution of silica Qn units

    if Q4_k05[-1] - p_Q4 - (p_B4*CQ4) < 0:
        next_Q4 = 0
    else:
        next_Q4 = Q4_k05[-1] - p_Q4 - (CQ4*p_B4)

    if Q3_k05[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4) < 0:
        next_Q3 = 0
    else:
        next_Q3 = Q3_k05[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4)

    if Q2_k05[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3) < 0:
        next_Q2 = 0
    else:
        next_Q2 = Q2_k05[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3)

    if Q1_k05[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2) < 0:
        next_Q1 = 0
    else:
        next_Q1 = Q1_k05[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2)

    if Q0_k05[-1] + p_Q1 + (p_B4*CQ1) < 0:
        next_Q0 = 0
    else:
        next_Q0 = Q0_k05[-1] + p_Q1 + (p_B4*CQ1)

    B3_k05.append(next_B3)
    B4_k05.append(next_B4)
    B2_k05.append(next_B2)
    B1_k05.append(next_B1)
    B0_k05.append(next_B0)

    Q4_k05.append(next_Q4)
    Q3_k05.append(next_Q3)
    Q2_k05.append(next_Q2)
    Q1_k05.append(next_Q1)
    Q0_k05.append(next_Q0)
    
    
plt.plot(M2O_k05, Q4_k05, 'r-', M2O_k05, Q3_k05, 'b-', M2O_k05, Q2_k05, 'g-', M2O_k05, Q1_k05, 'y-', M2O_k05, Q0_k05, 'k-')
plt.plot(mod_si_data_k05, Q4_data_k05, 'rd', mod_si_data_k05, Q3_data_k05, 'bd', mod_si_data_k05, Q2_data_k05, 'gd', mod_si_data_k05, Q1_data_k05, 'yd', mod_si_data_k05, Q0_data_k05, 'kd')
plt.axis([0, 65, 0, 65])
plt.xlabel("Modifier mol %")
plt.ylabel(f"Qn species concentration")
plt.title('Qn distribution K = 4')
plt.savefig(os.path.join(f'{fil}_data', f'{fil} Silicate_k4_0w.png'))
plt.show()

plt.plot(
#        M2O_k05, B3_k05, 'r-', 
        M2O_k05, B4_k05, 'b-'
#        , M2O_k05, B2_k05, 'g-', M2O_k05, B1_k05, 'y-', M2O_k05, B0_k05, 'k-'
        )
plt.plot(mod_b_data_k05, N4_data_k05, 'bd')
plt.axis([0, 65, 0, 65])
plt.xlabel("Modifier mol %")
plt.ylabel(f"Bn species concentration")
plt.title('Bn distribution K = 4')
plt.savefig(os.path.join(f'{fil}_data', f'{fil} Borate_k4_0w.png'))
plt.show()

m_Qmod_k4 = np.column_stack([M2O_k05, Q4_k05, Q3_k05, Q2_k05, Q1_k05, Q0_k05])
np.savetxt(os.path.join(f'{fil}_data', f"k4_Qn__model_{fil}.csv"), m_Qmod_k4)

m_Qdata_k4 = np.column_stack([mod_si_data_k05, Q4_data_k05, Q3_data_k05, Q2_data_k05, Q1_data_k05, Q0_data_k05])
np.savetxt(os.path.join(f'{fil}_data', f"k4_Qn__data_{fil}.csv"), m_Qdata_k4)

m_Bmod_k4 = np.column_stack([M2O_k05, B4_k05])
np.savetxt(os.path.join(f'{fil}_data', f"k4_Bn_model_{fil}.csv"), m_Bmod_k4)

m_Bdata_k4 = np.column_stack([mod_b_data_k05, N4_data_k05])
np.savetxt(os.path.join(f'{fil}_data', f"k4_Bn_data_{fil}.csv"), m_Bdata_k4)



Q4_mp_k6 = []
Q3_mp_k6 = []
Q2_mp_k6 = []
Q1_mp_k6 = []
Q0_mp_k6 = []

N4_mp_k6 = []

k_number = 6
draw_range = (((1/(k_number+1))*300)+((k_number/(1+k_number))*400))
draw_nr = list(range(int(draw_range)))
draw_ar = np.array(draw_nr)

M2O_k6 = []

for i2 in draw_ar:
    next_mod = (((draw_ar[i2]) / ((k_number*0.5/(1+(k_number*0.5)))*200 + ((1 / (1 + (k_number*0.5)))*100) + draw_ar[i2]))*100)
    M2O_k6.append(next_mod)

M2O_k6.append(75)

B4_B2 = []
for i3 in draw_ar:
    if M2O_k6[i3] < res.x[1]:
        B4_B2.append(1)
    else:
        B4_B2.append(0)

B4_B2 = np.array(B4_B2)

B3_k6 = [(1/((k_number*0.5)+1))*100, ]
B4_k6 = [0, ]
B2_k6 = [0, ]
B1_k6 = [0, ]
B0_k6 = [0, ]


Q4_k6 = [((k_number*0.5)/(1+(k_number*0.5)))*100, ]
Q3_k6 = [0, ]
Q2_k6 = [0, ]
Q1_k6 = [0, ]
Q0_k6 = [0, ]

# The function that makes each iteration at a time


for i4 in draw_ar:


    p_Q4 = Q4_k6[-1]*w[4] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])
    p_Q3 = Q3_k6[-1]*w[5] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])
    p_Q2 = Q2_k6[-1]*w[6] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])
    p_Q1 = Q1_k6[-1]*w[7] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])
    
    p_B3 = B3_k6[-1]*w[0] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])
    p_B4 = B4_k6[-1]*w[1] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])
    p_B2 = B2_k6[-1]*w[2] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])
    p_B1 = B1_k6[-1]*w[3] / (Q4_k6[-1]*w[4] + Q3_k6[-1]*w[5] + Q2_k6[-1]*w[6] + Q1_k6[-1]*w[7] + B3_k6[-1]*w[0] + B4_k6[-1]*w[1] + B2_k6[-1]*w[2] + B1_k6[-1]*w[3])

    #Contribution to N4 from Si
    CQ4 = p_Q4 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ3 = p_Q3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ2 = p_Q2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CQ1 = p_Q1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

    #Contribution to N4 from B
    CB3 = p_B3 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CB2 = p_B2 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)
    CB1 = p_B1 / (p_Q4 + p_Q3 + p_Q2 + p_Q1 + p_B3 + p_B2 + p_B1)

    # Evolution of borate Qn units
    if B3_k6[-1] - p_B3 - (p_B4 * CB3) < 0:
        next_B3 = 0
    else:
        next_B3 = B3_k6[-1] - p_B3 - (p_B4 * CB3)

    if B4_k6[-1] + (p_B3*B4_B2[i4]) - p_B4 < 0:
        next_B4 = 0
    else:
        next_B4 = B4_k6[-1] + (p_B3*B4_B2[i4]) - p_B4

    if B2_k6[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2) < 0:
        next_B2 = 0
    else:
        next_B2 = B2_k6[-1] + (p_B3*(1-B4_B2[i4])) + (p_B4 * CB3) + p_B4 - p_B2 - (p_B4 * CB2)

    if B1_k6[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
        next_B1 = 0
    else:
        next_B1 = B1_k6[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

    if B0_k6[-1] + p_B1 + (p_B4 * CB1 ) < 0:
        next_B0 = 0
    else:
        next_B0 = B0_k6[-1] + p_B1 + (p_B4 * CB1)


    #Now for the evolution of silica Qn units

    if Q4_k6[-1] - p_Q4 - (p_B4*CQ4) < 0:
        next_Q4 = 0
    else:
        next_Q4 = Q4_k6[-1] - p_Q4 - (CQ4*p_B4)

    if Q3_k6[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4) < 0:
        next_Q3 = 0
    else:
        next_Q3 = Q3_k6[-1] - p_Q3 + p_Q4 + (p_B4*CQ4) - (CQ3*p_B4)

    if Q2_k6[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3) < 0:
        next_Q2 = 0
    else:
        next_Q2 = Q2_k6[-1] - p_Q2 + p_Q3 - (p_B4*CQ2) + (p_B4*CQ3)

    if Q1_k6[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2) < 0:
        next_Q1 = 0
    else:
        next_Q1 = Q1_k6[-1] - p_Q1 + p_Q2 - (p_B4*CQ1) + (p_B4*CQ2)

    if Q0_k6[-1] + p_Q1 + (p_B4*CQ1) < 0:
        next_Q0 = 0
    else:
        next_Q0 = Q0_k6[-1] + p_Q1 + (p_B4*CQ1)

    B3_k6.append(next_B3)
    B4_k6.append(next_B4)
    B2_k6.append(next_B2)
    B1_k6.append(next_B1)
    B0_k6.append(next_B0)

    Q4_k6.append(next_Q4)
    Q3_k6.append(next_Q3)
    Q2_k6.append(next_Q2)
    Q1_k6.append(next_Q1)
    Q0_k6.append(next_Q0)
    
    
plt.plot(M2O_k6, Q4_k6, 'r-', M2O_k6, Q3_k6, 'b-', M2O_k6, Q2_k6, 'g-', M2O_k6, Q1_k6, 'y-', M2O_k6, Q0_k6, 'k-')
plt.plot(mod_si_data_k6, Q4_data_k6, 'rd', mod_si_data_k6, Q3_data_k6, 'bd', mod_si_data_k6, Q2_data_k6, 'gd', mod_si_data_k6, Q1_data_k6, 'yd', mod_si_data_k6, Q0_data_k6, 'kd')
plt.axis([0, 65, 0, 65])
plt.xlabel("Modifier mol %")
plt.ylabel(f"Qn species concentration")
plt.title('Qn distribution K = 6')
plt.savefig(os.path.join(f'{fil}_data', f'{fil} Silicate_k6_0w.png'))
plt.show()

plt.plot(
#        M2O_k05, B3_k05, 'r-', 
        M2O_k6, B4_k6, 'b-'
#        , M2O_k05, B2_k05, 'g-', M2O_k05, B1_k05, 'y-', M2O_k05, B0_k05, 'k-'
        )
plt.plot(mod_b_data_k6, N4_data_k6, 'bd')
plt.axis([0, 65, 0, 65])
plt.xlabel("Modifier mol %")
plt.ylabel(f"Bn species concentration")
plt.title('Bn distribution K = 6')
plt.savefig(os.path.join(f'{fil}_data', f'{fil} Borate_k6_0w.png'))
plt.show()

m_Qmod_k6 = np.column_stack([M2O_k6, Q4_k6, Q3_k6, Q2_k6, Q1_k6, Q0_k6])
np.savetxt(os.path.join(f'{fil}_data', f"k6_Qn__model_{fil}.csv"), m_Qmod_k6)

m_Qdata_k6 = np.column_stack([mod_si_data_k6, Q4_data_k6, Q3_data_k6, Q2_data_k6, Q1_data_k6, Q0_data_k6])
np.savetxt(os.path.join(f'{fil}_data', f"k6_Qn__data_{fil}.csv"), m_Qdata_k6)

m_Bmod_k6 = np.column_stack([M2O_k6, B4_k6])
np.savetxt(os.path.join(f'{fil}_data', f"k6_Bn_model_{fil}.csv"), m_Bmod_k6)

m_Bdata_k6 = np.column_stack([mod_b_data_k6, N4_data_k6])
np.savetxt(os.path.join(f'{fil}_data', f"k6_Bn_data_{fil}.csv"), m_Bdata_k4)
