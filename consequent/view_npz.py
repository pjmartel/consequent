import numpy as np
from random import sample, seed
#from matplotlib  import use
# use('Agg')
import matplotlib.pyplot as plt
from sys import argv
from scipy.stats import gumbel_r, norm


var_dict = np.load(argv[1])['saved_vars'].tolist()
sscores = var_dict['sscores']
uscore = var_dict['uscore']
N = var_dict['N']


# Fit extreme value distribution to data
miu, beta = gumbel_r.fit(sscores)
m, v = norm.fit(sscores)

print("Length of sscores: ", len(sscores))
print("Calculed histogram for {} scramble scores".format(N))
print("Max scrambled score:", max(sscores))
print("Min scrambled score:", min(sscores))
print("Unscrambled score:", uscore)
print("Median of scrambled scores:", np.median(sscores))
print("Gumbel miu:", miu)
print("Gumbel beta:", beta)
print()
# print("Aligment matrix:")
# np.savetxt(sys.stdout, ma, fmt="%3d")

# basename = "smith_{}_{}_{}_{:3.1f}_{:3.1f}".format(
#    N, len(seqA), matrix, abs(gapOpen), abs(gapExtend))
fig, ax = plt.subplots(figsize=(17, 9))
if(uscore):
    # ax.set_xticks(list(np.arange(0,160,10))+[uscore])
    ax.set_xticks(list(np.arange(10, 70, 10)))
    ax.tick_params(labelsize=15)
    # plt.gca().get_xticklabels()[-1].set_color('red')
    # ax.set_xlim(0,100)
# ax.set_title("S-W, {} aligns,len {}, matrix {}, gapo {}, gape {}".format(
#    N, len(seqA), matrix, gapOpen, gapExtend))
counts, bins, _ = ax.hist(sscores, bins=np.arange(
    min([10]), max(sscores), 1), align='left', rwidth=0.90)
x = np.arange(bins[0], bins[-1], 0.01)
ax.plot(x, sum(counts)*(bins[1]-bins[0])*gumbel_r.pdf(x, miu, beta), color="red", linewidth="3",
        label="miu: {:5.2f}, beta: {:5.2f}".format(miu, beta))
#plt.vlines(uscore, 0, max(counts),linewidth=8,color="green")
ax.legend()
plt.show()

#print("Saving plot to '"+basename+".pdf"+"'")
# plt.savefig(basename+".pdf")
#print("Saving data to '"+basename+".npy"+"'")
#np.save(basename, sscores)
