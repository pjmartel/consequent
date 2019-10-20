# align_check_mp.py --
#
# Aligns two sequences locally or globally and tests
# the signifcance of the aligment score bv comparing
# it with the distribution of scores arising from a set of
# random aligments. The random alignments are produced bv
# scrambling one of the sequences repeteadly.
# The final distribution of random scores is fitted to
# an extreme value (Gumbel) distribution, and the
# fitted CDF is used to calculate the probability of obtaining
# the aligment score in a random aligment -if that probability
# is very low, there is good indication that sequences are
# related, as it is very unlikely that two unrelated sequences
# would produce the same score.
#
# Note 1: the distribution of random scores is in fact not
# according to an extreme value distribution for a set of
# random global aligment scores. Moreover, the exact distribution
# is not known. As such, this method cannot be reliably used
# to estimate a signifcance value for a global alignment score.
#
# Note 2: scrambling one of the sequences as a means of simulating
# the aligment of two unrelated sequences is questionable on
# theoretical grounds, since sequences do not evolve and diverge
# by means of shufflings at the individual aminoacid level. It is
# expected the residual similarity between to distantly related
# sequences be higher than that of scrambled aligments. Therefore
# this method is likely overestimating the significance of the
# unscrambled aligment.
#
#                                                 ---PJM 2019
#
#  Pre-requisites: numpy, sciypy, matplotlib, Biopython, tqdm
#
#  This is the multiprocessing enabled version.
#
import numpy as np
from random import sample, seed
# from matplotlib import use
# use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gumbel_r
from matrix import readScoreMatrix, getMatrix
from alignment import smithWaterman, needlemanWunsch
from sequence import getUniprotSeq
from multiprocessing import Process, Manager
import argparse
from tqdm import tqdm


def scrambler_aligner(pn, ssd, N, sa, sb, ms, go, ge):
    seed()  # Without this, the random sequences are the same
    # for every process!
    sscores = []
    # Show the tqdm progress bar only for process 0 !
    for i in tqdm(range(N)) if pn == 0 else range(N):
        # print("Process {}, pass {} ".format(pn,i+1)
        sb = "".join(sample(sb, len(sb)))  # scramble b sequence
        s, a, ma, ta = alignFunction(
            sa, sb, ms, gapO=go, gapE=ge, ScoreOnly=True)

        sscores.append(s)

    ssd[pn] = sscores


# seqB = "HEAGAWGHEE"
# seqA = "PAWHEAE"
# seqB = "GVTAH"
# seqA = "AVTLI"
seqA = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
seqB = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
# seqB = "MVHYKLMCFDVRGLGEVIRQLFYLGDVSFEDFRVSREEFKSLKSNLPSGQLPVLEIDGVM"
gapOpen = -10.0
gapExtend = -1.0
smatrix = "BLOSUM62"


# Parse commmand line arguments
parser = argparse.ArgumentParser(
    "smithwat_mp.py", description="Sequence aligment significance analysis.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-A", "--alignment-method",
                    help="Alignment method (local or global).", choices=['global', 'local'],
                    default="local", type=str)
parser.add_argument("-N", "--random-scrambles",
                    help="Number of random scrambles.", default=100, type=int)
parser.add_argument("-o", "--gap-opening",
                    help="Penalty score for gap opening.", default=-10, type=float)
parser.add_argument("-e", "--gap-extension",
                    help="Penalty socre for gap extension.", default=-1, type=float)
parser.add_argument("-S", "--sequences", nargs=2, metavar=('UniprotA', 'UniprotB'),
                    help="Sequences to align (Uniprot Codes)",
                    type=str, default=["", ""])
parser.add_argument("-m", "--score-matrix",
                    help="Scoring matrix file.", default="BLOSUM62", type=str)
parser.add_argument("-p", "--cores",
                    help="Number of cores to use in computation.", default=4, type=int)
parser.add_argument("-s", "--seed",
                    help="Random number generator seed.", default=1234, type=int)


print()
print("/// Statistical analysis of sequence alignments. ///")
print()
print()

args = parser.parse_args()
N = args.random_scrambles
core_count = args.cores
gapOpen = -args.gap_opening
gapExtend = -args.gap_extension
smatrix = args.score_matrix
seed(args.seed)
if args.sequences[0] != "":
    seqA = getUniprotSeq(args.sequences[0])
    print("--> ", seqA.description)
    seqB = getUniprotSeq(args.sequences[1])
    print("--> ", seqB.description)
    seqA = str(seqA.seq)
    seqB = str(seqB.seq)

print()
if(args.alignment_method == "local"):
    print("Using Smith-Waterman local alignment algorithm.")
    alignFunction = smithWaterman
else:
    print("Using Needleman-Wunsh global alignment algorithm.")
    alignFunction = needlemanWunsch

# print(N, gapOpen, gapExtend, smatrix)
# print(seqA)
# print(seqB)


readScoreMatrix(smatrix)

matScore = getMatrix()

# Calculate unscrambled alignment and score
s, a, ma, ta = alignFunction(
    seqA, seqB, matScore, gapO=gapOpen, gapE=gapExtend)
ua = a
uscore = s

print("Scoring matrix: ", smatrix)
print("Unscrambled score:", uscore)
print("Unscrambled identity: {:.2%}".format(sum([ua[0][i] == ua[1][i] and
                                                 ua[0][i] != '-' for i in range(len(ua[0]))])/len(ua[0])))
print("Alignment length:", len(ua[0]))
print("Unscrambled alignment:")
print()
w = 60  # alignment width
for i in range(1 + len(ua[0]) // w):
    print("SeqA - ", ua[0][w*i:w*(i+1)], w*(i+1) % len(ua[0]))
    print("SeqB - ", ua[1][w*i:w*(i+1)], w*(i+1) % len(ua[0]))
    print()
print()

if N == 0:
    exit(0)

print("Calculating distribution of scores for {} scrambled alignments.".format(N))
print("Using {} processor cores.".format(core_count))

# Distribute random scrambles across cores
N_per_core = [N//core_count + (1 if k < (N % core_count) else 0)
              for k in range(core_count)]
#core_count = 4
print("Steps per core: ", N_per_core)

# Setup up multiprocessing
procs = []
sscores_dict = Manager().dict()

# Launch subprocesses
for i in range(core_count):
    proc = Process(target=scrambler_aligner, args=(
        i, sscores_dict, N_per_core[i], seqA, seqB, matScore, gapOpen, gapExtend))
    procs.append(proc)
    proc.start()

for proc in procs:
    proc.join()

# print(sscores_dict.values())
sscores = sum(sscores_dict.values(), [])


# Fit extreme value distribution to the scramble alignment data
miu, beta = gumbel_r.fit(sscores)
print("Length of sscores: ", len(sscores))
print("Computing histogram for {} scramble scores".format(N))
print("Max scrambled score:", max(sscores))
print("Min scrambled score:", min(sscores))
print("Median of scrambled scores:", np.median(sscores))
print("Gumbel miu:", miu)
print("Gumbel beta:", beta)
print("Probability of unscrambled score in a random alignment: ",
      1-gumbel_r.cdf(uscore, miu, beta))
print()

# Generate the basename for save files
basename = "smith" if args.alignment_method == "local" else "needle"
basename += "_{}_{}_{}_{:3.1f}_{:3.1f}".format(
    N, len(seqA), smatrix, abs(gapOpen), abs(gapExtend))

# Create the plot
fig, ax = plt.subplots()
ax.set_title("S-W, {} aligns,len {}, matrix {}, gapo {}, gape {}".format(
    N, len(seqA), smatrix, gapOpen, gapExtend))
counts, bins, _ = ax.hist(sscores, bins=np.arange(
    min(sscores), max(sscores)), align='left', rwidth=0.95)
x = np.arange(bins[0], bins[-1], 0.01)
ax.plot(x, sum(counts)*(bins[1]-bins[0])*gumbel_r.pdf(x, miu, beta),
        color="red", linewidth="3",
        label="miu: {:5.2f}, beta: {:5.2f}".format(miu, beta))
ax.legend()

# Save the plot to a pdf file
print("Saving plot to '"+basename+".pdf"+"'")
plt.savefig(basename+".pdf")

# Save the data to an npz compressed dictionary

# Extract a subdictionary from locals()
#  containing only the variable we need to save
L = locals()  # can't use "locals()" inside dict comprehension
saved_vars = {k: L[k] for k in ('N',
                                'gapOpen',
                                'gapExtend',
                                'smatrix',
                                'uscore',
                                'sscores',
                                'seqA',
                                'seqB')}

# Add some extra variables
saved_vars['alignment_function'] = alignFunction.__name__
saved_vars['_description'] = ""  # add a comment field

# Save it to a compressed npz file
print("Saving data to '"+basename+".npz"+"'")
np.savez_compressed(basename, saved_vars=saved_vars)
