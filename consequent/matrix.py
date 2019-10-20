import numpy as np
from sys import exit

mat = np.zeros((26,26), dtype = np.int8)   


def OLD_readScoreMatrix(filename):
    try:
        with open(filename, "r") as f:
            f = open(filename, "r")
            t = f.read().splitlines()
    except IOError:
        print("Can't open scoring matrix file '{}'.".format(filename))
        exit(1)

    t[-1] = t[-1].strip().split()
    global mat
    for i in range(len(t)-1):
        t[i] = t[i].split()
        # print(t[i])
        # print(t[i][0]," ",end="")
        for j in range(1, len(t[i])):
            # print(t[i][j]," ",end="")
            mat[ord(t[i][0])-65, ord(t[-1][j-1])-65] = t[i][j]
        # print()
    mat[:] = mat + mat.transpose() - np.diag(np.diag(mat))


def readScoreMatrix(filename):

    try:
        with open(filename) as f:
            lines = (line for line in f if not line.startswith('#'))
            M = np.loadtxt(lines, dtype=int, skiprows=1, converters ={0: lambda s : ord(s)-ord('A')})  
    except FileNotFoundError:
        print("Score Matrix file '{}' not found. Exiting.".format(filename))
        exit(1)
    #M = np.loadtxt("BLOSUM50",dtype=int, skiprows=7,
    #                converters ={0: lambda s : ord(s)-65})
    M = M[:-1,:-1]
    for i in range(26): 
        for j in range(26): 
            if i in M[:,0] and j in M[:,0]: 
                mat[i,j] = M[M[:,0] == i][0][1+np.where(M[:,0] == j)[0][0]]


def sPairScore(a, b):
    p1 = ord(a) - ord('A')
    p2 = ord(b) - ord('A')
    if p1 > 26 or p1 < 0 or p2 > 26 or p2 < 0 :
        print("Error: invalid aminoacid code in sPairScore.")
        exit(1)
    
    return mat[p1, p2]


def nPairScore(a,b):
    return mat[a, b]


def getMatrix():
    return mat

