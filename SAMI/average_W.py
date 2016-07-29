# simple script to calculate average velocity width of a.70 catalogue

import numpy as np

alfalfa = open('a70_160624.csv')
alfalfa_lines = alfalfa.readlines()

for i in range(len(alfalfa_lines)):
    alfalfa_lines[i] = alfalfa_lines[i].split(',')

alfalfa.close()

W = []
for i in range(1,len(alfalfa_lines)):
    W.append(float(alfalfa_lines[i][7]))

print len(W)
print np.mean(W)
