import scipy
import numpy as np
from scipy import spatial
from geopy.distance import great_circle

sami = open('SAMI_EarlyDataRelease.txt')
sami_lines = sami.readlines()

for i in range(len(sami_lines)):
    sami_lines[i] = sami_lines[i].split(' ')

sami.close()

alfalfa = open('a70_160624.csv')
alfalfa_lines = alfalfa.readlines()

for i in range(len(alfalfa_lines)):
    alfalfa_lines[i] = alfalfa_lines[i].split(',')

alfalfa.close()

sami_coord = []

for i in range(len(sami_lines)):
    if i != 0 and i != 108:
        try:
            sami_coord.append((float(sami_lines[i][3]),float(sami_lines[i][7])))
        except:
            sami_coord.append((float(sami_lines[i][3]),float(sami_lines[i][8])))

alfalfa_coord = []

for i in range(len(alfalfa_lines)):
    if i != 0:
        alfalfa_coord.append((float(alfalfa_lines[i][4]),float(alfalfa_lines[i][5])))

closest = scipy.spatial.distance.cdist(sami_coord, alfalfa_coord)

#for i in range(len(closest)):
#    print min(closest[i]), sami_coord[i], alfalfa_coord[np.where(closest[i] == min(closest[i]))[0]]

## Other method

r = 1

new_closest = []

for i in range(len(sami_coord)):
    new_closest.append([])
    for j in range(len(alfalfa_coord)):
        theta1 = sami_coord[i][0]*np.pi/180.0
        phi1 = sami_coord[i][1]*np.pi/180.0
        theta2 = alfalfa_coord[j][0]*np.pi/180.0
        phi2 = alfalfa_coord[j][1]*np.pi/180.0
        x1 = np.sin(phi1)*np.cos(theta1)
        y1 = np.sin(phi1)*np.sin(theta1)
        z1 = np.cos(phi1)
        x2 = np.sin(phi2)*np.cos(theta2)
        y2 = np.sin(phi2)*np.sin(theta2)
        z2 = np.cos(phi2)
#        dist = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
#                dist = ((np.sin(sami_coord[i][1]*np.pi/180.0)*np.cos(sami_coord[i][0]*np.pi/180.0) - np.sin(alfalfa_coord[j][1]*np.pi/180.0)*np.cos(alfalfa_coord[j][0]*np.pi/180.0))**2+(np.sin(sami_coord[i][1]*np.pi/180.0)*np.sin(sami_coord[i][0]*np.pi/180.0) - np.sin(alfalfa_coord[j][1]*np.pi/180.0)*np.sin(alfalfa_coord[j][0]*np.pi/180.0))**2+(np.cos(sami_coord[i][1]*np.pi/180.0) - np.cos(alfalfa_coord[j][1]*np.pi/180.0))**2)**0.5
    #    beta = np.arccos(1 - dist**2/2)*180.0/np.pi
        #beta = np.arccos(np.cos(phi1)*np.cos(phi2)+np.sin(phi1)*np.sin(phi2)*np.cos(theta1 - theta2))*180/np.pi
###        beta = np.arctan2((np.cos(phi2)*(np.sin(theta1-theta2))**2+(np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*(np.cos(theta1-theta2))**2))**0.5,np.sin(phi1)*np.sin(phi2)+np.cos(phi1)*np.cos(phi2)*np.cos(theta1-theta2))
        beta = great_circle((theta1,phi1),(theta2,phi2)).meters*180.0/(np.pi*6371000)*3600
        new_closest[i].append(beta)

new_closest = np.array(new_closest)

#print new_closest
for i in range(len(new_closest)):
    print min(new_closest[i]), sami_coord[i], alfalfa_coord[np.where(new_closest[i] == min(new_closest[i]))[0]]
