"""
Script to produce protein and nucleic acid dot plots.

-- PJMartel 2018

"""
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, Circle
import numpy as np
from sys import argv
from score_matrix import readScoreMatrix, getMatrix
from get_seq import getUniprotSeq
import matplotlib.pyplot as plt
import matplotlib
from scipy.signal import convolve2d as cv2
from matplotlib.ticker import MultipleLocator
from matplotlib.widgets import Slider
import argparse

class DotPlot:
    def __init__(self, score_matrix='blosum50.txt', seqs=('ADE','WCA')):
        """
        Set the background and grid.
        Initialize the parent Figure class.
        """

        readScoreMatrix(score_matrix)
        print("Using scoring matrix '{}'".format(score_matrix))
        self.seqs = seqs

    def createPlot(self, window_size=1, cut_off=-1):
        # Convert seqs to numerical
        seqA = np.array(list(map(ord, self.seqs[0])), dtype=np.int8)
        seqB = np.array(list(map(ord, self.seqs[1])), dtype=np.int8)
        seqA -= ord("A")
        seqB -= ord("A")
        plot_size = (len(seqA), len(seqB))
        (ix, iy) = np.meshgrid(range(len(seqB)), range(len(seqA)))
        self.dotMatrix = getMatrix()[seqB[ix], seqA[iy]]
        M = self.dotMatrix

        print("Min and Max befor kernel: ", M.min(), M.max())

        if window_size > 1 :
            kernel = np.eye(window_size)
            M = cv2(M, kernel, mode='same', boundary='fill')
            trim = (window_size-1)//2
            M = M[trim:-trim, trim:-trim]
            print("Min and Max after kernel: ", M.min(), M.max())

        if cut_off != -1 :
            c = M.min() + (M.max()-M.min())*cut_off
            #print("Min, Max, c-par, cut_off: ", M.min(), M.max(), c, cut_off)
            #print(len(M[M <= c]))
            M[ M <= c ] = 0
            M[ M != 0 ] = 1
            #print("Max and min elements;", M.min(), M.max())
            #print("Count of 0's and 1's :", M.size-M.sum(), M.sum(), (M.size-M.sum())+M.sum() )

        return M






    
                      
    # def calcScore(self, posA, posB, seqs, win_size):
    #     if posA < win_size or posA > len(seqs[0])-win_size:
    #         return 0
    #     if posB < win_size or posB > len(seqs[1])-win_size:
    #         return 0
    #     score = 0
    #     for i in range(-win_size, +win_size):
    #         score += nPairScore(seqs[0][i], seqs[1][i])
    #     return score


    
def main():

    def update(val, cut_off, window_size):
        if cut_off != -1 :
            cut_off = scut.val
        #else :
        #    custom_cmap = cmap.
        print("Cutoff ", scut.val)
        window_size = int(swin.val)
        M = A.createPlot(window_size, cut_off)
        im.set_data(M)
        #im.set_cmap('')

    # Parse commmand line arguments
    parser = argparse.ArgumentParser(
        "dot_plot.py", description="Biological sequence dotplot comparator.")
    parser.add_argument("-w", "--window-size",
                        help="Size of the averaging window.", default=1, type=int)
    parser.add_argument("-c", "--cut-off",
                        help="Cutoff for plot coloring.", default=-1, type=float)
    parser.add_argument("-s", "--size", nargs=2, metavar=('width', 'height'),
                        help="Plot size (width height)",
                        type=int, default=[10, 10])
    parser.add_argument("-m", "--score-matrix",
                        help="Scoring matrix file.", default="blosum50.txt", type=str)

    parser.add_argument('UniprotA', help="Uniprot sequence code")
    parser.add_argument('UniprotB', help="Uniprot sequence code")

    args = parser.parse_args()
    width = args.size[0]
    height = args.size[1]
    
    window_size = args.window_size
    cut_off = args.cut_off
    plot_size = args.size
    score_matrix = args.score_matrix
    
    seqA = getUniprotSeq(args.UniprotA)
    seqB = getUniprotSeq(args.UniprotB)

    print("Window size: ", window_size)
    print("Cutoff: ", cut_off)
    print("Plot size: ", (width, height))   
    print(seqA.description)
    print(seqB.description)

    A = DotPlot(score_matrix=score_matrix, seqs=(seqA,seqB))
    M = A.createPlot(window_size, cut_off)

    #plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    #plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

    matplotlib.rc('axes',labelsize=15)
    #fig, ax = plt.subplots(figsize=(10.6,10.6))
    fig, ax = plt.subplots(figsize=(10.6,10.6))
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator(50))
    #ax.set_xticks(major_ticks)
    #ax.set_xticks(minor_ticks, minor=True)
    #ax.set_yticks(major_ticks)
    #ax.set_yticks(minor_ticks, minor=True)
    
    ## disable tick labels
    #labels = [item.get_text() for item in ax.get_xticklabels()]
    #empty_string_labels = ['']*len(labels)
    #ax.set_xticklabels(empty_string_labels)
    #labels = [item.get_text() for item in ax.get_yticklabels()]
    #empty_string_labels = ['']*len(labels)
    #ax.set_yticklabels(empty_string_labels)
    
    ## enable grid
    #ax.grid(which='both')
    #ax.grid(which='minor', alpha=0.2)
    #ax.grid(which='major', alpha=0.5)

    ## set labels   
    #ax.set_xlabel(seqA.description)
    #ax.set_ylabel(seqB.description)
    ax.set_xlabel(args.UniprotA)
    ax.set_ylabel(args.UniprotB)
    #ax.set_xlabel(seqA.seq)
    #ax.set_ylabel(seqB.seq[::-1])

    ## move the x axis indicators to top
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position('top')

    ## sliders
    axcolor = 'lightgoldenrodyellow'
    axcut = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    axwin = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    scut = Slider(axcut, 'Cut-off', 0.0, 1.0, valinit=cut_off)
    swin = Slider(axwin, 'Window-size', 1, 25, valinit=window_size, valstep=2)
    scut.on_changed(lambda e: update(e, cut_off, window_size))
    swin.on_changed(lambda e: update(e, cut_off, window_size))

    #if cut_off == -1 :
    #  axcmap = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    #  scmap = Slider(axcmpa, 'Gradient', 0, 0.4, valinit=window_size)
    ## show it
    im = ax.imshow(M, cmap='Greys', interpolation='none')
    fig.colorbar(im, ax=ax, shrink=0.8)
    plt.show()
    
if __name__ == "__main__" :
    main()
