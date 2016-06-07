#!/usr/bin/python

###
# an example of how to extract information from txt tables
###

file_1 = open('KMOS_GCLASS_CLUS0034.txt', 'r')

lines = file_1.readlines()

for i in range(len(lines)):
    lines[i] = lines[i].split(' ')
    for j in range(len(lines[i])):
        if j < 3:
            lines[i][j] = int(lines[i][j])
        else:
            lines[i][j] = float(lines[i][j])

file_1.close()

print lines
