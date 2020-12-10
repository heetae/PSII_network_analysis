# Copyright (C) 2020 by
# Heetae Kim <heetaekim@udd.cl>
# All rights reserved.
# GNU General Public License v3.0.

import pandas as pd
import numpy as np

def Gini(list_x):
    x=sortlist(list_x)[0]
    fx=sortlist(list_x)[1]
    return GiniIndex(trapezoide(x,fx))

def sortlist(list_x): # turn a list to x coordinates and values
    sort_list=np.sort(list_x)
    cummulative_list=[]
    x=np.arange(0,len(list_x)+1)/len(list_x)
    temp_total=np.sum(list_x)
    cummulative_sum=0.0
    cummulative_list.append(cummulative_sum)
    if temp_total >0:
        for j in range(len(list_x)):
            cummulative_sum+=sort_list[j]
            cummulative_list.append(cummulative_sum/float(temp_total))
    else: # temp_total=0
        for j in range(len(list_x)):
            cummulative_list.append(0)
    return x,cummulative_list

def trapezoide(x, fx):
    out=0
    if np.count_nonzero(fx) == 0:
        out=0.5
    if not np.count_nonzero(fx) ==1:
        for i in range(len(x)-1):
            h = (x[i+1] - x[i])
            out += h * (fx[i+1] + fx[i]) / 2
    return(out)

def GiniIndex(x): # get the Gini index
    A=0.5-x
    gini= A/0.5
    return gini

def shannon(ary):
    log2ary=np.log2(ary)
    return -1*np.sum(ary*log2ary)

def analyze(FRET_RATE_FILE):
    # data load
    data=pd.read_csv('{}'.format(FRET_RATE_FILE))

    samples=[str(i[0]) for i in np.array(data.iloc[:,:1])]
    idx2char = sorted(list(set(samples)))  # index -> char
    idx2char.append('Photosynthesis') # virtual node for photosynthesis, id 314
    idx2char.append('Intrinsic dissipation') # virtual node for self dissipation, id 315
    char2idx = {c: i for i, c in enumerate(idx2char)}

    Q=np.zeros((len(idx2char),len(idx2char)), dtype=float)
    # Construct the matrix Q
    for i in data.index:
        row=str(data.iloc[i,0])
        column=str(data.iloc[i,1])
        FRETrate=float(data.iloc[i,2])
        Q[char2idx[row]][char2idx[column]]=FRETrate
    centers={}
    for j,i in enumerate(idx2char):
        if i in ['A_405','D_402','a_405','d_402']:
            centers[i]=j

    # Additional rates
    for i in [centers['a_405'],centers['d_402'],centers['A_405'],centers['D_402']]:
        Q[i][len(Q)-2]=1./1.5 # charge separation rate
    for i in range(len(Q)-2): # except for the virtual nodes
        Q[i][len(Q)-1]=0.0005 # Intrinsic dissipation

    # Fill the diagonal elements
    N=len(Q)
    for i in range(N):
        OutTransient=0
        OutTransient+=sum(Q[i,:i]) # from 0 to i-1
        OutTransient+=sum(Q[i,i+1:]) # from i+1 to end
        Q[i,i]=-OutTransient    

    # Convert the transient matrix to the probability matrix, P    
    I=np.identity(len(idx2char), dtype=float)
    l=min(np.diag(Q))
    P=I-Q/l         

    # Iteration
    A=np.ones((1,N)) # Initial excitation to chlorophylls
    A[0][-2]=0
    A[0][-1]=0

    energy_t={} # time series energy concentration of each molecule. index is time
    temp1=0
    temp2=1

    Gini_n=[]
    Gini_n_nophoto=[]
    for t in range(100000):
        energy_t[t]=[]
        for i in range(N):
            energy_t[t].append(A[0][i])
        Gini_n.append(Gini(energy_t[t]))
        Gini_n_nophoto.append(Gini(energy_t[t][:-2])) # only between chlorophylls
        A=np.matmul(A,P)

    energy_chlo={} # index is chl.
    for i in range(len(energy_t[0])):
        energy_chlo[i]=[]
        for j in range(len(energy_t)):
            energy_chlo[i].append(energy_t[j][i])
        
    DFenergy_t=pd.DataFrame()
    chlo=len(energy_t[0])
    for c in range(len(energy_chlo)):
        DFenergy_t[c]=energy_chlo[c]
    DFenergy_t['Gini_all']=Gini_n   # Gini coefficient including the virtual nodes
    DFenergy_t['Gini_chl']=Gini_n_nophoto    # Gini coefficient between chlorophylls, used in the study

    DFenergy_t.to_pickle('DFenergy_t_{}.pkl'.format(FRET_RATE_FILE[:-4]))


DataFiles = ['FRET_all_a.csv','FRET_major_a.csv','FRET_minor_a.csv','FRET_natural.csv','FRET_minor_b.csv','FRET_major_b.csv','FRET_all_b.csv']
for FRET_RATE_FILE in DataFiles:
    analyze(FRET_RATE_FILE)
