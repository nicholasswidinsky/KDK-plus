import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path

def ReadInFileNaI(file, numChannel):
    
    data = []
    for c in range(numChannel):
        data.append([ [] for _ in range(3)])
    # channels = [2*n for n in range(numChannel)]
    channels = [0,2,4,8]
    
    
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            for c in range(numChannel):
                if int(lines[1]) == channels[c]:
                    data[c][0].append(float(lines[2]))
                    data[c][1].append(float(lines[3]))
                    data[c][2].append(lines[5])
    return data

def sumChannels(file):
    
    i = 0
    sumInt  = 0 
    sumTotals = []
    times = []
    timeDiffs = []
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            sumInt+=int(lines[1])
            times.append(int(lines[2])/1000)
            
            i+=1
            
            if i == 4:
                sumTotals.append(sumInt)
                timeDiffs.append(times[-1]-times[0])
                i = 0
                sumInt = 0
                times = []
            
    return sumTotals, timeDiffs
                
            
            

filepath = "/home/nick/PhD/KDK+/Quad_coinc_testing/two_func_generator/not_synchronized/"

GammaDetector = 'NaI'
file = "SDataR_quad_coinc_test_2_func_500mV_3kHz_5kHz_50ns_short_coinc_window_corrected.csv"
# file = 'SDataR_Large_LSC_vessel_NaI_Co60_quad_coinc_Vertical_scatter_v4.CSV'


timeData = []
with open(f'{filepath}/{file}') as f:
    next(f)
    for line in f:
        lines = line.split('\n')[0].split(';')
        
        timeData.append(float(lines[2]))
        
        

fileData = ReadInFileNaI(f'{filepath}/{file}', 4)
sumTotals,timeDiffs = sumChannels(f'{filepath}/{file}')


print(len(fileData[0][0]))
print(len(fileData[1][0]))
print(len(fileData[2][0]))
print(len(fileData[3][0]))

# print(timeDiffs[0:100])
# print(f'Sum of groups of 4 channels:')


for i in range(len(sumTotals)):
    if sumTotals[i] != 14:
        print(4*i)
        
        print((timeData[4*i]-timeData[0])*10**(-12))
        break