import numpy as np
from random import sample, seed
#import matplotlib.pyplot as plt
from sys import argv, stdout
#from scipy.stats import gumbel_r
from score_matrix import readScoreMatrix, getMatrix
from seqali import smithWaterman, smithFast, plotMat, plotTraceMat
from multiprocessing import Process, Manager

def scrambler_aligner(pn, ssd, N, sa, sb, ms, go, ge):
    seed()
    sscores = []
    for i in range(N):
        #print("Process {}, pass {} ".format(pn,i+1))
        sa = "".join(sample(sa, len(sa)))
        s, a, ma, ta = smithFast(
            sa, sb, ms, gapO=go, gapE=ge)

        sscores.append(s)
        
    ssd[pn] = sscores 
        

#seqB = "HEAGAWGHEE"
#seqA = "PAWHEAE"
# seqB = "GVTAH"
# seqA = "AVTLI"
seqB = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
seqA = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
#seqB = "MVLSPADKTNVKAAWGKVGAHAGEYG"
#seqA = "MVHLTPEEKSAVTALWGKVNVDEVGG"

gapOpen = -10
gapExtend = -1
#gapOpen = -8
#gapExtend = -8
matrix = "BLOSUM50"

if(len(argv) > 1):
    N = int(argv[1])
else:
    N = 100

    
# init score matrix
#matScore = np.zeros((26, 26), dtype=np.int8)
#readMat("blosum50.txt", matScore)

readScoreMatrix(matrix)

matScore = getMatrix()

# Calculate unscrambled aligment and score
s, a, ma, ta = smithWaterman(
    seqA, seqB, matScore, gapO=gapOpen, gapE=gapExtend)
ua = a
uscore = s

print("Scoring matrix: ", matrix)
print("Unscrambled score:", uscore)
print("Unscrambled identity: {:.2%}".format(sum([ua[0][i] == ua[1][i] and
                               ua[0][i] != '-' for i in range(len(ua[0]))])/len(ua[0])))
print("Unscrambled alignment:")
print("SeqA - ", ua[0])
print("SeqB - ", ua[1])
print()

if N==0 :
    exit(0)
    
print("Calculating distribution of scrambled alignment scores.")


proc_count = 4
procs = []
sscores_dict = Manager().dict()

for i in range(proc_count):
    proc = Process(target=scrambler_aligner, args=(i, sscores_dict, N, seqA, seqB, matScore, gapOpen, gapExtend))
    procs.append(proc)
    proc.start()

for proc in procs:
    proc.join()

#print(sscores_dict.values())
sscores = sum(sscores_dict.values(),[])
#print(sscores)
#exit(0)
N = len(sscores) # for 4 cores its 4 times the initial value


# Fit extreme value distribution to data
#miu, beta = gumbel_r.fit(sscores)

print("Length of sscores: ", len(sscores))
print("Calculed histogram for {} scramble scores".format(N))
print("Max scrambled score:", max(sscores))
print("Min scrambled score:", min(sscores))
print("Median of scrambled scores:", np.median(sscores))
print("Gumbel miu:", miu)
print("Gumbel beta:", beta)
print()
# print("Aligment matrix:")
# np.savetxt(sys.stdout, ma, fmt="%3d")

print("Saving data to","'smith_{}_{}_{}_{:3.1f}_{:3.1f}.npy'".format(
    N, len(seqA), matrix, abs(gapOpen), abs(gapExtend)))
np.save("smith_{}_{}_{}_{:3.1f}_{:3.1f}".format(
    N, len(seqA), matrix, abs(gapOpen), abs(gapExtend)),sscores)
