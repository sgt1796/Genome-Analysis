import scipy as sp
import numpy as np
import sys
from matplotlib import pyplot as plt

class FastaAna:
    def __init__(self, filename):
        self.filename = filename
        self.seq = self.openFASTA(filename)
        self.length = len(self.seq)

    def openFASTA(self, filename):
        file = open(filename,"r")
        line = file.readline()
        text = file.read().replace("\n","")
        if line[0] == ">":
            return text
        else:
            text += line[:-1]
        return text

    def statSeq(self):
        seq = self.seq
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

    def findSTART(self, start, end):
        # find all start codons and return position
        seq = self.seq
        count = 0
        pos = []
        #r = input("enter the range of search (if thorough enter 0,-1):")
        #start = int(r.split(",")[0])
        #end = int(r.split(",")[1])
        if(end == -1):
            end = self.length - 1
        for x in range(start,end):
            if(seq[x] == 'A' and seq[x+1] == 'T' and seq[x+2] == 'G'):
                pos.append(x)
                count += 1
                x += 2
        #print("total " + str(count)+ " START found:")
        #print(pos)
        return pos

    def findSTOP(self, start, end):
        #find all stop codons and return position
        seq = self.seq
        count = 0
        pos = []
        #r = input("enter the range of search (if thorough enter 0,-1):")
        #start = int(r.split(",")[0])
        #end = int(r.split(",")[1])
        if(end == -1):
            end = self.length - 1
        for x in range(start,end):
            stop = False
            if(seq[x] == 'T'):
                if(seq[x+1] == 'A'):
                    if(seq[x+2] == 'G' or seq[x+2] == 'A'):
                        stop = True
                if(seq[x+1] == 'G'):
                    if(seq[x+2] == 'A'):
                        stop = True
            if(stop):
                pos.append(x)
                count += 1
                x += 2
        #print("total " + str(count)+ " STOP found:")
        #print(pos)
        return pos

    def findORF(self, fr):
        #not fix: need to add frame
        pStart = []
        pStop = []
        for x in self.findSTART(0,-1):
            if(x%3 == fr-1):
                pStart.append(x)
        for x in self.findSTOP(0,-1):
            if(x%3 == fr-1):
                pStop.append(x)
        #print(pStart)
        #print(pStop)
        hold, count = 0, 0
        ORF = []
        for x in pStart:
            if(x > pStop[hold]):
                for y in range(hold, len(pStop)-1):
                    if(x < pStop[y]):
                        ORF.append((x,pStop[y]))
                        hold = y
                        count += 1
                        break

        print("Total "+str(count)+" ORF found in reading frame " + str(fr) + ":")
        print(ORF)
        return ORF

    def plotORF(self, *ORF):
        x1 = np.arange(self.length)
        y1 = x1
        for orf in ORF:
            x,y = [],[]
            for w in orf:
                x.append(w[0])
                y.append(w[1])
                plt.plot([w[0],w[0]], [w[0],w[1]])
            plt.plot(x,y,'.',x1,y1,'--')
        plt.show()








covFile = r"/Users/GuotaiShen/Desktop/Bioinformatic/2019-nCoV/sequence.fasta"
file1 =  r"/Users/GuotaiShen/Desktop/Bioinformatic/pNMT1-GNAS/Pnmt1-GNAS.txt"
CoV = FastaAna(covFile)
orf1 = CoV.findORF(1)
orf2 = CoV.findORF(2)
orf3 = CoV.findORF(3)
CoV.plotORF(orf1, orf2, orf3)
