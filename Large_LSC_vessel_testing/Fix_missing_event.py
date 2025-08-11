import numpy as np
import matplotlib.pyplot as plt


def ReadInFileNaI(file, numChannel):
    
    data = []
    for c in range(numChannel):
        data.append([ [] for _ in range(4)])
    # channels = [2*n for n in range(numChannel)]
    channels = [0,2,6,8]
    
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            for c in range(numChannel):
                if int(lines[1]) == channels[c]:
                    data[c][0].append(int(lines[1]))
                    data[c][1].append(float(lines[2]))
                    data[c][2].append(float(lines[3]))
                    data[c][3].append(lines[5])
    return data

def fixMissingData(data):
    j = 0
    
    
    for i in range(len(data[0][0])):
        sum = data[0][0][i] + data[1][0][i] + data[2][0][i] + data[3][0][i]
        print(sum)
        if sum != 16:
            print(i)
            
        j += 1
        if j == 4:
            j = 0
        


filepath = '/home/nick/PhD/KDK+/Large_LSC_testing/Quad_coinc/2025_03_21/'

# GammaDetector = 'LYSO'
# file = 'SDataR_Small_LSC_vessel_LYSO_Cs137_triple_coinc_Vertical_scatter_v4_10CG_LSC.CSV'
# bckfile = '/home/nick/PhD/KDK+/Small_LSC_testing/2025_03_19/SDataR_Small_LSC_vessel_LYSO_bck_triple_coinc_Vertical_scatter_v4_10CG_LSC.CSV'
# # bckfile = f'{filepath}/{file}'

# fileData = ReadInFile(f'{filepath}/{file}', 4)
# filebck = ReadInFile(f'{bckfile}', 4)

GammaDetector = 'NaI'
file = 'SDataR_Large_LSC_vessel_NaI_Co60_quad_coinc_Vertical_scatter_v4.CSV'
bckfile = '/home/nick/PhD/KDK+/Small_LSC_testing/2025_03_19/SDataR_Small_LSC_vessel_NaI_bck_triple_coinc_Vertical_scatter_v4_10CG_LSC.CSV'
        
fileData = ReadInFileNaI(f'{filepath}/{file}', 4)
filebck = ReadInFileNaI(f'{bckfile}', 4)

fixMissingData(fileData)