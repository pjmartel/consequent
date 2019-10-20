import numpy as np
import matplotlib.pyplot as plt
from sys import exit


def plotMat(mat, fs):
    fig, axs = plt.subplots()
    axs.axis("tight")
    axs.axis("off")
    the_table = axs.table(cellText=mat, cellLoc='center', loc='center')
    for key, cell in the_table.get_celld().items():
        row, col = key
        cell.set_fontsize(fs)
    fig.show()


def plotTraceMat(tmat, smat, fs):
    fig, axs = plt.subplots()
    axs.axis("tight")
    axs.axis("off")
    arrows = [8598, 8593, 8592]
    plotMat = [['' for c in range(tmat.shape[1])]
               for r in range(tmat.shape[0])]
    for i in range(tmat.shape[0]):
        for j in range(tmat.shape[1]):
            # if(smat[i, j] == 0):
            #    plotMat[i][j] = ' '
            # else:
            if smat[i, j] != 0:
                plotMat[i][j] = chr(arrows[tmat[i, j]])
    # plotMat = np.fromfunction(lambda i, j: arrows[mat[i, j]], mat.shape)
    t1 = axs.table(cellText=plotMat, cellLoc='center', loc='center')
    t1.set_fontsize(fs)
    for key, cell in t1.get_celld().items():
        cell.set_linewidth(0)

    # t2 = axs.table(cellText=smat, cellLoc='center', loc='center')
    # t2.set_fontsize(fs)
    # t2.set_alpha(0.5)

    fig.show()


def needlemanWunsch(seqA, seqB, scoreMat, gapO=-10, gapE=-0.5, ScoreOnly=False):
    """
Computes the optimal global aligment of two amino acid
sequences using the Needleman-Wunsch algorithm.

seqA, seqB - aminoacid sequences
scoreMat - scoring matrix
gapO - gap opening penalty
gapE - gap extension penalty
ScoreOnly - if True, only the optimal score is computed
            (aligment is returned as "None")
    """

    # Check for valid score matrix argument
    if not isinstance(scoreMat, np.ndarray):
        print("Error in 'smitWaterman': scoring matrix has wrong data type.")
        exit(0)
    elif scoreMat.shape != (26, 26):
        print("Error in 'smitWaterman': scoring matrix has wrong shape.")
        exit(0)

    # print(type(scoreMat)
    # Convert seqs to numerical
    seqA = np.array(list(map(ord, seqA)), dtype=np.int8)
    seqB = np.array(list(map(ord, seqB)), dtype=np.int8)
    seqA -= ord("A")
    seqB -= ord("A")

    matAli = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int32)
    matTrace = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int8)

    # init aligment and trace matrices
    matAli[1:, 0] = np.fromfunction(
        lambda i: gapO+gapE*i, (len(seqA),))
    matTrace[1:, 0] = 1
    #
    # (ix,iy) = np.meshgrid(range(len(seqB)), range(len(seqA)))
    # matAli[1:
    matAli[0, 1:] = np.fromfunction(
        lambda j: gapO+gapE*j, (len(seqB),))
    matTrace[0, 1:] = 2
    # Let'a initialize the rest of the aligment matrix with
    # the comparison scores from the scoring matrix.
    # Note the meshgrid trick to avoid loops ;-)
    (ix, iy) = np.meshgrid(range(len(seqB)), range(len(seqA)))
    matAli[1:, 1:] = scoreMat[seqB[ix], seqA[iy]]

    # for i in range(len(A)+1):
    #     matAli[i, 0] = (gapOpen+gapExtend*(i-1))*(i > 0)
    #     matTrace[i, 0] = 1
    # for j in range(len(B)+1):
    #     matAli[0, j] = (gapOpen+gapExtend*(j-1))*(j > 0)
    #     matTrace[0, j] = 2
    # matTrace[0, 0] = 0

    # Fill aligment matrix
    for i in range(1, len(seqA)+1):
        for j in range(1, len(seqB)+1):
            # tempScores = (matAli[i-1, j-1] + matAli[i, j],
            tempScores = (matAli[i-1, j-1] + scoreMat[seqA[i-1], seqB[j-1]],
                          matAli[i-1, j] +
                          (gapE if (matTrace[i-1, j] == 1) else gapO),
                          matAli[i, j-1] + (gapE if (matTrace[i, j-1] == 2) else gapO))
            idx = np.argmax(tempScores)
            matTrace[i, j] = idx
            matAli[i, j] = tempScores[idx]

    Score = matAli[-1, -1]

    # if ScoreOnly is True, return before computing alignment
    if ScoreOnly:
        return(Score, None, matAli, matTrace)

    # Backtrack aligment
    alB = []
    alA = []
    i = len(seqA)
    j = len(seqB)
    cur = matTrace[i, j]
    while(i > 0 or j > 0):
        # print(i, j)
        cur = matTrace[i, j]
        if cur == 0:
            alA.append(chr(seqA[i-1]+65))
            alB.append(chr(seqB[j-1]+65))
            i -= 1
            j -= 1
        elif cur == 1:
            alA.append(chr(seqA[i-1]+65))
            alB.append("-")
            i -= 1
        elif cur == 2:
            alB.append(chr(seqB[j-1]+65))
            alA.append("-")
            j -= 1
        else:
            print("Invalid value in aligment matrix.")
            exit(0)

    Alignment = np.array(["".join(alA[::-1]), "".join(alB[::-1])])

    return(Score, Alignment, matAli, matTrace)


def smithWaterman(seqA, seqB, scoreMat, gapO=-10, gapE=-0.5, ScoreOnly=False):
    """
Computer the optimal local aligment of two aminoacid
sequences using the Smith-Waterman algoritm.

seqA, seqB - aminoacid sequences
scoreMat - scoring matrix
gapO - gap opening penalty
gapE - gap extension penalty
ScoreOnly - if True, only the optimal score is computed
            and no aligment is returned.
    """

    # Check for valid score matrix argument
    if not isinstance(scoreMat, np.ndarray):
        print("Error in 'smitWaterman': scoring matrix has wrong data type.")
        exit(1)
    elif scoreMat.shape != (26, 26):
        print("Error in 'smitWaterman': scoring matrix has wrong shape.")
        exit(1)

    # Convert seqs to numerical
    seqA = np.array(list(map(ord, seqA)), dtype=np.int8)
    seqB = np.array(list(map(ord, seqB)), dtype=np.int8)
    seqA -= ord("A")
    seqB -= ord("A")

    matAli = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int32)
    matTrace = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int8)
    # matTrace = np.zeros((len(A)+1,len(B)+1),dtype=(np.int8,2))
    # matTrace = np.empty((len(A)+1,len(B)+1),dtype=object)

    # init aligment and trace matrices
    # matAli[1:, 0] = np.fromfunction(
    #    lambda i: gapO+gapE*i, (len(seqA),))
    # matTrace[1:, 0] = 1
    # matAli[0, 1:] = np.fromfunction(
    #    lambda j: gapO+gapE*j, (len(seqB),))
    # matTrace[0, 1:] = 2

    # for i in range(len(A)+1):
    #     matAli[i, 0] = (gapOpen+gapExtend*(i-1))*(i > 0)
    #     matTrace[i, 0] = 1
    # for j in range(len(B)+1):
    #     matAli[0, j] = (gapOpen+gapExtend*(j-1))*(j > 0)
    #     matTrace[0, j] = 2
    # matTrace[0, 0] = 0

    # Fill aligment matrix
    for i in range(1, len(seqA)+1):
        for j in range(1, len(seqB)+1):
            tempScores = (matAli[i-1, j-1] + scoreMat[seqA[i-1], seqB[j-1]],
                          matAli[i-1, j] +
                          (gapE if (matTrace[i-1, j] == 1) else gapO),
                          matAli[i, j-1] + (gapE if (matTrace[i, j-1] == 2) else gapO), 0)
            idx = np.argmax(tempScores)
            matAli[i, j] = tempScores[idx]
            if idx == 3:
                idx = 2
            matTrace[i, j] = idx

    Score = matAli.max()

    if ScoreOnly:
        return(Score, None, matAli, matTrace)

    # Backtrack aligment
    alB = []
    alA = []
    # (i, j) are the indices of max score in matAli
    i, j = np.unravel_index(np.argmax(matAli), matAli.shape)

    # print("--->", Score, i, j)
    # print(matAli)
    while(matAli[i, j] != 0):
        # print(i, j)
        cur = matTrace[i, j]
        if cur == 0:
            alA.append(chr(seqA[i-1]+65))
            alB.append(chr(seqB[j-1]+65))
            i -= 1
            j -= 1
        elif cur == 1:
            alA.append(chr(seqA[i-1]+65))
            alB.append("-")
            i -= 1
        elif cur == 2:
            alB.append(chr(seqB[j-1]+65))
            alA.append("-")
            j -= 1
        else:
            print("Invalid value in trace matrix:", cur)
            exit(0)

    Alignment = np.array(["".join(alA[::-1]), "".join(alB[::-1])])

    return(Score, Alignment, matAli, matTrace)


def needleFast(seqA, seqB, scoreMat, gapO=-8, gapE=-8):

    # Convert seqs to numerical
    seqA = np.array(list(map(ord, seqA)), dtype=np.int8)
    seqB = np.array(list(map(ord, seqB)), dtype=np.int8)
    seqA -= ord("A")
    seqB -= ord("A")

    matAli = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int32)
    matTrace = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int8)

    # init aligment and trace matrices
    matAli[1:, 0] = np.fromfunction(
        lambda i: gapO+gapE*i, (len(seqA),))
    matTrace[1:, 0] = 1
    #
    # (ix,iy) = np.meshgrid(range(len(seqB)), range(len(seqA)))
    # matAli[1:
    matAli[0, 1:] = np.fromfunction(
        lambda j: gapO+gapE*j, (len(seqB),))
    matTrace[0, 1:] = 2
    # Let'a initialize the rest of the aligment matrix with
    # the comparison scores from the scoring matrix.
    # Note the meshgrid trick to avoid loops ;-)
    (ix, iy) = np.meshgrid(range(len(seqB)), range(len(seqA)))
    matAli[1:, 1:] = scoreMat[seqB[ix], seqA[iy]]

    # Fill aligment matrix
    for i in range(1, len(seqA)+1):
        for j in range(1, len(seqB)+1):
            # tempScores = (matAli[i-1, j-1] + matAli[i, j],
            tempScores = (matAli[i-1, j-1] + scoreMat[seqA[i-1], seqB[j-1]],
                          matAli[i-1, j] +
                          (gapE if (matTrace[i-1, j] == 1) else gapO),
                          matAli[i, j-1] + (gapE if (matTrace[i, j-1] == 2) else gapO))
            idx = np.argmax(tempScores)
            matTrace[i, j] = idx
            matAli[i, j] = tempScores[idx]
    return(matAli[-1, -1], None, matAli, matTrace)


def smithFast(seqA, seqB, scoreMat, gapO=-8, gapE=-8):

    # Convert seqs to numerical
    seqA = np.array(list(map(ord, seqA)), dtype=np.int8)
    seqB = np.array(list(map(ord, seqB)), dtype=np.int8)
    seqA -= ord("A")
    seqB -= ord("A")

    matAli = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int32)
    matTrace = np.zeros((len(seqA)+1, len(seqB)+1), dtype=np.int8)

    # Fill aligment matrix
    for i in range(1, len(seqA)+1):
        for j in range(1, len(seqB)+1):
            tempScores = (matAli[i-1, j-1] + scoreMat[seqA[i-1], seqB[j-1]],
                          matAli[i-1, j] +
                          (gapE if (matTrace[i-1, j] == 1) else gapO),
                          matAli[i, j-1] + (gapE if (matTrace[i, j-1] == 2) else gapO), 0)
            idx = np.argmax(tempScores)
            matAli[i, j] = tempScores[idx]
            if idx == 3:
                idx = 2
            matTrace[i, j] = idx

    Score = matAli.max()

    return(Score, None, matAli, matTrace)
