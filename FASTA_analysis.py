import scipy as sp
import numpy as np
import sys
from matplotlib import pyplot as plt

class FastaAna:
    def __init__(self, filename):
        self.filename = filename
        self.seq = self.openFASTA(filename)
        self.revS = self.revSeq(self.seq)
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

    def findSTART(self, start, end, rev = False):
        # find all start codons and return position
        seq = self.seq
        if (rev):
            seq = self.revS
        count = 0
        pos = []
        #r = input("enter the range of search (if thorough enter 0,-1):")
        #start = int(r.split(",")[0])
        #end = int(r.split(",")[1])
        if(end == -1):
            end = self.length - 1
        while start < end-2:
            if(seq[start] == 'A' and seq[start+1] == 'T' and seq[start+2] == 'G'):
                pos.append(start)
                count += 1
                start += 2
            start += 1
        #print("total " + str(count)+ " START found:")
        #print(pos)
        return pos

    def findSTOP(self, start, end, rev = False):
        #find all stop codons and return position
        seq = self.seq
        if(rev):
            seq = self.revS
        count = 0
        pos = []
        #r = input("enter the range of search (if thorough enter 0,-1):")
        #start = int(r.split(",")[0])
        #end = int(r.split(",")[1])
        if(end == -1):
            end = self.length - 1
        while start < end-2:
            stop = False
            if(seq[start] == 'T'):
                if(seq[start+1] == 'A'):
                    if(seq[start+2] == 'G' or seq[start+2] == 'A'):
                        stop = True
                if(seq[start+1] == 'G'):
                    if(seq[start+2] == 'A'):
                        stop = True
            if(stop):
                pos.append(start)
                count += 1
                start += 2
            start += 1
        #print("total " + str(count)+ " STOP found:")
        #print(pos)
        return pos

    def findORF(self, fr = 1):
        #There are 6 frames in total: 1, 2, 3, -1, -2, -3
        rev = False
        if fr < 0:
            rev = True
            fr = -fr

        #find all start and stop condons and store them in pStart and pStop
        pStart = []
        pStop = []
        for x in self.findSTART(0,-1,rev):
            if(x%3 == fr-1):
                pStart.append(x)
        for x in self.findSTOP(0,-1,rev):
            if(x%3 == fr-1):
                pStop.append(x)
        #print(pStart)
        #print(pStop)

        #find ORF by the shortest distance
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

        if(rev):
            print("Total " + str(count) + " ORF found in reading frame -" + str(fr) + ":")
        else:
            print("Total " + str(count) + " ORF found in reading frame " + str(fr) + ":")
        print(ORF)
        return ORF

    def plotORF(self, threshold, *ORF, coor = 2000):
        print(ORF)
        if(not ORF):
            ORF = [self.findORF(1)]
        x1 = np.arange(self.length)
        y1 = x1
        count = 1
        for orf in ORF:
            x,y = [],[]
            for w in orf:
                if(w[1]-w[0] > threshold):
                    x.append(w[0])
                    y.append(w[1])
                    plt.plot([w[0],w[0]], [w[0],w[1]])
                    if w[1]-w[0] > coor:
                        plt.text(w[0]+10, w[1]+10, ' [%d-%d]\n %d nt'%(w[0], w[1], w[1]-w[0]))

            plt.plot(x,y,'.' ,label = 'reading frame %d'%count)
            plt.plot(x1, y1, '--')
            plt.legend(loc='lower right')
            count += 1
        plt.show()

    def revSeq(self, seq):
        #reverse seq
        rev = ""
        for x in seq:
            if(x == 'A'):
                rev += 'T'
            elif(x == 'T'):
                rev += 'A'
            elif(x == 'C'):
                rev += 'G'
            else:
                rev += 'C'
        ans = rev[::-1]
        return ans


    def getSeq(self, start, end):
        return self.seq[start:end]






#covFile = r"/Users/GuotaiShen/Desktop/Bioinformatic/2019-nCoV/sequence.fasta
covFile = r'C:\Users\Guotai Shen\Desktop\Genenome seq\2019-nCoV.fasta'
#file1 =  r"/Users/GuotaiShen/Desktop/Bioinformatic/pNMT1-GNAS/Pnmt1-GNAS.txt"
CoV = FastaAna(covFile)
orf1 = CoV.findORF(1)
orf2 = CoV.findORF(2)
orf3 = CoV.findORF(3)
Rorf1 = CoV.findORF(-1)
Rorf2 = CoV.findORF(-2)
Rorf3 = CoV.findORF(-3)
#CoV.plotORF(100, orf1, orf2, orf3, Rorf1, Rorf2, Rorf3)
