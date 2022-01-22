import numpy as np
import matplotlib.pyplot as plt
import pydot
import pickle
import pathlib
import os
from scipy import stats
import random

def single_dim_lineage(namein,dictin,tier = False):
    gong = ['AB','C','D','MS','E']
    Subsetylist = [[],[],[],[],[]]
    Subsetxlist = [[],[],[],[],[]]
    namelist,yaxis,xaxis = sorted(list(dictin), key = len ),[],[];target_xaxis=[];target_yaxis=[]
    if True:
        for n in namelist:
            for m in range(0,len(gong)):
                if gong[m] in n:
                    Subsetylist[m].append(dictin[n])
                    Subsetxlist[m].append(len(re.sub(r'[A-Z]', '', n)))
                    continue
    namein = namein.split('/')[-1][:-4]
    for n in range(0,len(Subsetxlist)):
        plt.scatter(np.array(Subsetxlist[n]),Subsetylist[n], label=gong[n], s = [30] * len(Subsetxlist[n]), alpha = .35)
    plt.title('Cellular Division Timing Events in WT C. Elegans Embryo\n'+namein)
    plt.xlabel('Division Event Count');plt.ylabel('Division Event Timing Increments')
    plt.legend(loc='upper left');plt.show()
    
def dict_to_scatter_lineage(name,lis,x,y,img = False,seps='', tier = False):
    gong = ['AB','C','D','MS','E']
    Subsetylist = [[],[],[],[],[]]
    Subsetxlist = [[],[],[],[],[]]
    xaxis = [];yaxis = []
    target_xaxis=[];target_yaxis=[]
    for n in (list(set(lis[x]).intersection(set(lis[y])))):
        for m in range(0,len(gong)):
            if gong[m] in n:
                Subsetylist[m].append(lis[x][n]);Subsetxlist[m].append(lis[y][n])
                continue
        xaxis.append(lis[x][n]);yaxis.append(lis[y][n])
    target_xaxis = np.array(target_xaxis);xaxis = np.array(xaxis)
    m, b, r, p_value, std_err = stats.linregress(np.append(xaxis,target_xaxis),yaxis+target_yaxis)
    #switch to siegelslope & theilsslope functions?
    name1,name2 = name[x].split('/')[-1][:-4],name[y].split('/')[-1][:-4]
    plt.title('Cellular Division Timing Discrepancies in WT C. Elegans Embryos')
    plt.xlabel(name1+'\nDivision Event Timing Increments');plt.ylabel(name2+'\nDivision Event Timing Increments')
    for n in range(len(Subsetxlist)):
        plt.scatter(np.array(Subsetxlist[n]),Subsetylist[n], label=gong[n], s = [30] * len(Subsetxlist[n]), alpha = .35)
    plt.legend(loc='upper left')
    plt.show()
    
def aggregate_rank_function(x):
    #reflist = sorted(x[0][0])
    reflist = sorted(x[0].tolist()[0]);ranklist = []
    for n in reflist:
        gimp = 0
        for m in x:
            k = m.tolist()[0]; gimp = gimp + k.index(n)
        ranklist.append((n,gimp))
    ranklist.sort(key=lambda x:x[1])
    return ranklist

def dist_seqences(l1,l2):
    try: 
        l1 = l1.tolist(); l2 = l2.tolist()
        if len(l1) == 1: l1 = l1[0]
        if len(l2) == 1: l2 = l2[0]
    except:
        if len(l1) == 1: l1 = l1[0]
        if len(l2) == 1: l2 = l2[0]
        #pass
    ref = []
    if len(l1) != len(l2): return
    for n in range(0,len(l1)): ref.append( (n,abs(l1.index(n)-l2.index(n)), int( (l1.index(n)+l2.index(n)) * 0.5 ) ))
    ref.sort(key=lambda x:(x[1],x[2]) )
    return ref

#hierarchical clustering on distance map, see if they align with date/time/occasion of imaging.
#build global time point data (add time with parents node time). Start from latest coordinated time point
#color code dots according to lineage molynobus? distance


    
