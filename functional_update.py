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
    
'''
    
WT_full_tree = [[ 3, 25,  4,  7,  9, 23, 22, 26, 17, 16, 24, 15, 29,  6,  2, 11,
          12, 21, 19, 13, 14,  5, 27, 20, 10, 18,  8, 28,  1,  0],
 [25, 23,  4,  3,  7, 16, 24, 22,  9, 17,  6, 18,  1, 19, 29, 26,
           2, 10, 28, 12, 15, 11,  0, 21,  5,  8, 20, 13, 14, 27],
 [ 3,  4, 17, 25, 23,  7, 26, 24, 19,  9, 28, 16, 11, 20,  0,  8,
           2, 10, 22, 13,  1, 29, 21,  5, 18,  6, 12, 14, 27, 15],
 [25,  3,  4,  9, 22, 23,  7, 26, 29, 15, 24, 17,  6, 11, 21, 16,
           2, 12, 20,  1, 19, 10, 18,  5, 13, 14, 27,  8, 28,  0],
 [20,  3, 23, 25, 17,  2,  4,  1, 24, 22,  7, 28, 16, 19, 10,  8,
          11,  9, 14, 18, 26,  0,  6, 15,  5, 21, 12, 29, 13, 27]]

WT_no_dang = [[18,  9, 13,  3, 16,  8, 20, 29, 26, 17, 28, 25,  4, 24, 27, 11,
           2, 22,  6, 23,  5, 19, 14, 10,  1, 15,  7, 21,  0, 12],
 [23, 25, 22, 18,  5, 17,  1, 21,  2, 27, 20, 16, 13,  6,  9, 24,
          15, 19, 26,  4, 10,  3,  7, 11, 29, 12, 28,  8,  0, 14],
 [ 3,  4, 25,  7, 23, 17, 11, 19, 24, 22,  2, 13, 16, 26, 10, 29,
           8, 21, 28,  1, 14, 27, 12,  5,  0,  9, 15, 18, 20,  6],
 [20,  1, 25, 17, 18, 19,  8, 13, 29, 26,  9, 10,  3,  5, 11, 27,
          16,  2,  4, 28, 15, 24, 14, 23, 12,  6,  0, 22, 21,  7],
 [20,  2, 17,  3,  1, 23, 19, 25,  4,  5, 24, 26,  9, 10,  8, 16,
          18, 27, 22, 28, 11, 21, 29, 14, 13, 15,  6,  0, 12,  7]]

WT_EXP_full_tree = [[ 3,  4,  9, 17, 26, 25, 19,  7, 20, 21, 18, 23,  8, 11,  0, 24,  1,
        12, 16, 15, 10,  5,  6, 28,  2, 14, 13, 29, 27, 22],
 [25, 23,  4,  3, 17, 11,  5, 22, 21, 20, 13, 12, 27, 15, 18, 29,  9,
         0,  8,  6, 26, 19,  1, 14, 28,  2, 24, 10, 16,  7],
 [ 5, 16,  4,  9, 18, 19, 21,  3,  6, 26,  8, 27, 28, 24, 25, 22, 12,
         1,  2, 20,  7, 29,  0, 23, 13, 10, 11, 17, 14, 15],
 [25,  3,  4, 23, 11, 17, 18, 15,  9, 12, 26,  0, 19, 29,  1, 20,  5,
         6,  7, 21,  8, 13, 24, 10, 16,  2, 14, 28, 22, 27],
 [20, 17, 23,  1, 19, 22,  9,  8, 26, 10, 25, 16, 21, 15, 12,  6, 13,
        27, 29,  5, 18,  4,  7,  0, 11, 14, 28, 24,  3,  2]]

WT_EXP_no_dang = [[ 8,  3, 19, 25,  6, 24,  2,  0, 10,  4, 20,  1,  9, 17, 23, 26, 12,
        28, 21, 22, 18, 11, 27, 14,  7,  5, 13, 16, 15, 29],
 [23, 25, 17, 22, 27,  6, 19,  5, 11, 20, 13, 18,  4, 12,  0,  2, 24,
        29, 21, 28,  3,  1, 14,  9, 10, 26,  7,  8, 16, 15],
 [ 5, 16,  4,  9, 18, 21,  6,  3, 26, 24, 28, 19, 27, 25, 22,  1,  2,
        20,  8, 12,  0,  7, 29, 23, 13, 10, 11, 17, 14, 15],
 [20,  1, 17, 19, 10,  8, 26,  9, 16, 23,  6, 12,  7,  0, 21, 18, 22,
        15, 11, 24, 14, 28,  4, 13,  3, 25, 29,  2,  5, 27],
 [20, 17,  1, 23, 19, 26,  9, 16, 10,  8, 22, 15, 21, 12,  6, 13, 18,
         7, 29, 25,  0, 27, 11, 14, 28,  5, 24,  4,  3,  2]]



sorttypes = ['R^2 cycle time','R^2 birth time', 'T.E.D', 'B.E.D. Cycle Time','B.E.D. Birth Time']


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
    ref.sort(key=lambda x:x[1]+x[2] )
    for n in ref: print(n)
    #return ref

for n in range(0,len(WT_full_tree)):
  print('WT/WT ' + sorttypes[n])
  print('(embryo ID, Difference In position, average position)')
  dist_seqences(WT_full_tree[n],WT_no_dang[n])
  q = input()

for n in range(0,len(WT_EXP_full_tree)):
  print('WT/EXP ' + sorttypes[n])
  print('(embryo ID, Difference In position, average position)')
  dist_seqences(WT_EXP_full_tree[n],WT_EXP_no_dang[n])
  q = input()
  
'''
#Copy paste commented out code into IDE to run comparision

#hierarchical clustering on distance map, see if they align with date/time/occasion of imaging.
#build global time point data (add time with parents node time). Start from latest coordinated time point
#color code dots according to lineage molynobus? distance


    
