import numpy as np
import matplotlib.pyplot as plt

def readInFile(file):
    time, energy, flag = [],[],[]
    with open(file) as f:
        for i in range(1):
            next(f)
        for line in f:
            data = line.split(';')
            # print(data)
            time.append(float(data[2]))
            energy.append(float(data[3]))
            flag.append(data[5].split('\n')[0])
            
        time = np.asfarray(time, dtype = float)
        energy = np.asfarray(energy, dtype = float)
        flag = np.array(flag)


    return time,energy,flag


def readInFileCombined(file):
    ch0Time, ch2Time, ch0Energy, ch2Energy, ch0Flag, ch2Flag = [],[],[],[],[],[]
    
    with open(file) as f:
        next(f)
        for line in f:
            data = line.split('\n')[0].split(';')
            
            if int(data[1]) == 0:
                ch0Time.append(float(data[2]))
                ch0Energy.append(float(data[3]))
                ch0Flag.append(data[5])
                
            if int(data[1]) == 2:
                ch2Time.append(float(data[2]))
                ch2Energy.append(float(data[3]))
                ch2Flag.append(data[5])
            
    return ch0Time, ch2Time, ch0Energy, ch2Energy, ch0Flag, ch2Flag   
            

def sortTimeData(Ch0Data, Ch1Data):
    
    sortedData = [[],[],[]]
    timeSep = 100000 #10ns
    print(len(Ch0Data[0]))
    for i in range(len(Ch0Data[0])):
        if i < 10:
            minRange = 0
        else:
            minRange = i-10
        if i > len(Ch0Data[0]) - 10:
            maxRange = len(Ch1Data[0])-1
        else:
            maxRange = len(Ch1Data[0])-1#i+10
        # print(i)
        for j in range(minRange,maxRange):
            # print(j)
            # print(f'Index: {j}')
            # print(maxRange)
            timeDiff = np.absolute(Ch0Data[0][i] - Ch1Data[0][j])
            # print(f'Time Difference: {timeDiff}')
            if timeDiff <= timeSep:
                       
            # if Ch0Data[0][i] <= Ch1Data[0][j] + timeDiff and Ch0Data[0][i] >= Ch1Data[0][i] - timeDiff:
                sortedData[0].append(Ch0Data[0][i])
                sortedData[1].append(Ch0Data[1][i])
                sortedData[2].append(Ch0Data[2][i])
            #     break
            #     Ch0Data = [np.delete(Ch0Data[0],i),np.delete(Ch0Data[1],i),np.delete(Ch0Data[1],i)]
                
    return sortedData
            

            
            
########################################################################################################################################################################################################
####                                                                                                                                                                                                ####
####                                                                                                                                                                                                ####
####                                               Initial testing with no coincidence requirements. Didn't work that well in the end since                                                         ####
####                                               I had to manually find the coincidences.                                                                                                         ####
####                                                                                                                                                                                                ####
########################################################################################################################################################################################################

# filepath = '/home/nick/PhD/KDK+/Large_LSC_testing/Initial_testing_no_coinc/taped_vessel_1600Hv/'

# ch0Filename = 'Data_CH0@V1730S_27411_second_test_no_coincidence_1600HV_taped_vessel.CSV'
# ch1Filename = 'Data_CH1@V1730S_27411_second_test_no_coincidence_1600HV_taped_vessel.CSV'


# ch0Time, ch0Energy, ch0Flag = readInFile(f'{filepath}/{ch0Filename}')
# ch1Time, ch1Energy, ch1Flag = readInFile(f'{filepath}/{ch1Filename}')

# Ch0Data = sortTimeData([ch0Time,ch0Energy,ch0Flag],[ch1Time,ch1Energy,ch1Flag])
# print(len(Ch0Data[0]))
# Ch1Data = sortTimeData([ch1Time,ch1Energy,ch1Flag],Ch0Data)

# print('Data lengths')
# print(len(Ch0Data[1]))
# print(len(Ch1Data[1]))


# xRange = [0,600]

# bins = 500
# binRange = np.linspace(xRange[0], xRange[1],  bins)


# # ax[1,0].hist2d(Zndf.iloc[:,17],Zndf.iloc[:,8], bins = [np.arange(xRangePlas[0], xRangePlas[1], binWidth*5), np.arange(xRangeLYSO[0],xRangeLYSO[1], binWidth*5)], cmin = 1)
# # ax[1,1].hist2d(df.iloc[:,17],df.iloc[:,8], bins = [np.arange(xRangePlas[0], xRangePlas[1], binWidth*5), np.arange(xRangeLYSO[0],xRangeLYSO[1], binWidth*5)], cmin = 1)


# fig, ax = plt.subplots(1,1)

# histCountsZn, xbins, ybins, image = ax.hist2d(Ch0Data[1],Ch1Data[1], bins = [binRange,binRange], cmin = 0)


# ax.set_xlabel('Channel 0 Energy')
# ax.set_ylabel('Channel 1 Energy')
# plt.show()

# timeDiff = []
# for i in range(len(Ch0Data[0])):
#     timeDiff.append(Ch0Data[0][i] - Ch1Data[0][i])

# bins = 50
# timeRange = [-10,10]
# timeBins = np.linspace(timeRange[0],timeRange[1])

# fig,ax = plt.subplots(1,1)

# ax.hist(timeDiff, bins = timeBins)

# plt.show()

########################################################################################################################################################################################################
####                                                                                                                                                                                                ####
####                                                                                                                                                                                                ####
####                                               Second try using on-board coincidence.                                                                                                           ####
####                                                                                                                                                                                                ####
####                                                                                                                                                                                                ####
########################################################################################################################################################################################################



filepath = '/home/nick/PhD/KDK+/Large_LSC_testing/Initial_testing_coinc/Compass_coinc_no_free_write/'
filename = 'SData_second_test_no_coincidence_1600HV_taped_vessel_Coinc.CSV'


ch0Time, ch2Time, ch0Energy, ch2Energy, ch0Flag, ch2Flag = readInFileCombined(f'{filepath}/{filename}')

print(len(ch0Energy))
print(len(ch2Energy))

xRange = [0,600]

bins = 500
binRange = np.linspace(xRange[0], xRange[1],  bins)

fig, ax = plt.subplots(1,1)
ax.hist2d(ch0Energy,ch2Energy, bins = [binRange,binRange], cmin = 0)

ax.set_xlabel('Ch0 Energy Integral')
ax.set_ylabel('Ch2 Energy Integral')

plt.show()