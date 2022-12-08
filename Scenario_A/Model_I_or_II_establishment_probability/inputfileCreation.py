import os
import sys
import numpy as np


# Parameters

tEstV = np.linspace(0, 100, 100)    # for numerical results
# tEstV = np.linspace(0, 100, 13)     # for stochastic simulations
betaV = [1e-5]
# betaV = [1e-6, 1e-5, 1e-4]
# phiV = [0.5]
phiV = [0.1, 0.5, 0.9]
sMV = [0.1]
cV = [0]
sigmaV = [0]
s0V = [-0.1]
alphaV = [1e-6]
N0V = [1e6]



# Deleting all previous inputfiles
dir = 'inputfiles/'
for f in os.listdir(dir):
    os.remove(os.path.join(dir, f))



# Writing new inputfiles
count = 1
for beta in betaV:
    for phi in phiV:
        for sM in sMV:
            for c in cV:
                for sigma in sigmaV:
                    for s0 in s0V:
                        for alpha in alphaV:
                            for N0 in N0V:
                                for tEst in tEstV:
                                    filepath = dir + "input_" + str(count) + ".dat"
                                    with open(filepath, 'w') as f:
                                        f.write("%e #beta\n" % beta)
                                        f.write("%f #phi\n" % phi)
                                        f.write("%f #sM\n" % sM)
                                        f.write("%f #c\n" % c)
                                        f.write("%f #sigma\n" % sigma)
                                        f.write("%e #alpha\n" % alpha)
                                        f.write("%d #N0\n" % N0)
                                        f.write("%f #s0\n" % s0)
                                        f.write("%f #tEst\n" % tEst)
                                    count += 1



# Returning the number of files on cluster
print(count - 1)
