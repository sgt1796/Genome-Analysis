import scipy as sp
import numpy as np
import sys

def openFASTA(filename):
    file = open(filename,"r")
    line = file.readline()
    text = file.read().replace("\n","")
    if line[0] == ">":
        return text
    else:
        text += line[:-1]
    return text

def StatSeq(seq):
    A,T,C,G = 0,0,0,0
    count = 0
    for x in seq:
        if(x == 'A'):
            A += 1
        elif(x == 'T'):
            T += 1
        elif(x == 'C'):
            C += 1
        else:
            G += 1
        count += 1
    stat = "total nt: "+str(count)+"\n A: "+str(A)+"\n T: "+str(T)+"\n C: "+str(C)+"\n G: "+str(G)
    print(stat)


seq = openFASTA(r"C:\Users\Guotai Shen\Desktop\Genenome seq\2019-nCoV.fasta")
StatSeq(seq)