import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams["font.size"] = 15
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["mathtext.default"] = 'regular'
matplotlib.rcParams['lines.markersize'] = 3

def ReadInFile(file, numChannel):
    
    data = []
    for c in range(numChannel):
        data.append([ [] for _ in range(3)])
    channels = [2*n for n in range(numChannel)]
    
    
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
    
    
def BinHistograms(data, Range, numBins):
    binRange = np.linspace(Range[0],Range[1],numBins)
    hist,binedges = np.histogram(data,binRange)
    
    return hist, binedges
    
def TimeDifference(data):
    timeDiff0_2, timeDiff0_4, timeDiff2_4 =[],[],[]
    
    for i in range(len(data[0][0])):
        timeDiff0_2.append(data[0][0][i]-data[1][0][i])
        timeDiff0_4.append(data[0][0][i]-data[2][0][i])
        timeDiff2_4.append(data[1][0][i]-data[2][0][i])
        
    timeDiff0_2 = np.asfarray(timeDiff0_2,float)
    timeDiff0_4 = np.asfarray(timeDiff0_4,float)
    timeDiff2_4 = np.asfarray(timeDiff2_4,float)
        
    return timeDiff0_2,timeDiff0_4,timeDiff2_4


######### 
#########
# Initial testing with large coinc window
#########
#########


# filepath = '/home/nick/PhD/KDK+/Large_LSC_testing/Triple_coinc_testing_CAEN_DAQ/2025_01_23/'
# filename = 'SDataR_Test_tripple_coinc_1700HV_100LSB_2.CSV'
# bckfilename = 'SDataR_Test_tripple_coinc_1700HV_100LSB_bck.CSV'

######### 
#########
# Initial testing with small coinc window
#########
#########

filepath = '/home/nick/PhD/KDK+/Large_LSC_testing/Triple_coinc_testing_CAEN_DAQ/2025_01_27/'
filename = 'SDataR_Test_tripple_coinc_1700HV_750LSB_Cs137_closer_lead_shorter_coinc_time.CSV'
bckfilename = 'SDataR_Test_tripple_coinc_1700HV_750LSB_bck_closer_lead_shorter_coinc_time.CSV'



data = ReadInFile(f'{filepath}/{filename}',3)
Bckdata = ReadInFile(f'{filepath}/{bckfilename}',3)

bins = 100
binRange = [0,1000]

binsLYSO = 100
binRangeLYSO = [0,2000]


CsInt1, CsBins1 = BinHistograms(data[0][1],binRange,bins)
CsInt2, CsBins2 = BinHistograms(data[1][1],binRange,bins)
CsIntLYSO, CsBinsLYSO = BinHistograms(data[2][1],binRangeLYSO,binsLYSO)

bckInt1, bckBins1 = BinHistograms(Bckdata[0][1],binRange,bins)
bckInt2, bckBins2 = BinHistograms(Bckdata[1][1],binRange,bins)
bckIntLYSO, bckBinsLYSO = BinHistograms(Bckdata[2][1],binRangeLYSO,binsLYSO)


normalize = False
logScale = False

fig, ax = plt.subplots(1,3, figsize = (30,10))

ax[0].hist(CsBins1[:-1],bins =CsBins1,density = normalize,weights=CsInt1/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Cs-137 Ch 0 (LSC)', log = logScale)
ax[0].hist(bckBins1[:-1],bins =bckBins1,density = normalize,weights=bckInt1/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

ax[1].hist(CsBins2[:-1],bins =CsBins2,density = normalize,weights=CsInt2/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Cs-137 Ch 2 (LSC)', log = logScale)
ax[1].hist(bckBins2[:-1],bins =bckBins2,density = normalize,weights=bckInt2/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

ax[2].hist(CsBinsLYSO[:-1],bins =CsBinsLYSO,density = normalize,weights=CsIntLYSO/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Cs-137 Ch 4 (LYSO)', log = logScale)
ax[2].hist(bckBinsLYSO[:-1],bins =bckBinsLYSO,density = normalize,weights=bckIntLYSO/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

ax[0].set_title('Left PMT')
ax[1].set_title('Right PMT')
ax[2].set_title('LYSO')

ax[0].set_xlabel('Integral (ADC Channel)')
ax[1].set_xlabel('Integral (ADC Channel)')
ax[2].set_xlabel('Integral (ADC Channel)')

ax[0].set_ylabel('Counts/bin/s')
ax[1].set_ylabel('Counts/bin/s')
ax[2].set_ylabel('Counts/bin/s')

plt.savefig(f'{filepath}/figures/{filename.split('.csv')[0]}_1D_histogram.png', dpi = 400)
plt.show()


timeRange = [-20000,20000]
timeBins = (timeRange[1] - timeRange[0])//2000

LSCbinRange = np.linspace(binRange[0],binRange[1],bins)
LYSObinRange = np.linspace(binRangeLYSO[0],binRangeLYSO[1],binsLYSO)

timebinsRange = np.linspace(timeRange[0],timeRange[1],timeBins)

timeDiff0_2,timeDiff0_4,timeDiff2_4 = TimeDifference(data)

fig,ax = plt.subplots(3,3,figsize = (30,30))

ax[0,1].hist(timeDiff0_2, bins = timebinsRange)
ax[0,1].set_title('Time Difference Ch0 - Ch2')

ax[0,2].hist(timeDiff0_4, bins = timebinsRange)
ax[0,2].set_title('Time Difference Ch0 - Ch4')

ax[1,0].hist(-1*timeDiff0_2, bins = timebinsRange)
ax[1,0].set_title('Time Difference Ch2 - Ch0')

ax[1,2].hist(timeDiff2_4, bins = timebinsRange)
ax[1,2].set_title('Time Difference Ch2 - Ch4')

ax[2,0].hist(-1 *timeDiff0_4, bins = timebinsRange)
ax[2,0].set_title('Time Difference Ch4 - Ch0')

ax[2,1].hist(-1 * timeDiff2_4, bins = timebinsRange)
ax[2,1].set_title('Time Difference Ch4 - Ch2')


plt.show()

fig, ax = plt.subplots(1,4,figsize = (30,10))

ax[0].hist2d(data[0][1],timeDiff0_4, bins = (LSCbinRange,timebinsRange),cmin = 1)
ax[1].hist2d(data[1][1],timeDiff2_4, bins = (LSCbinRange,timebinsRange),cmin = 1)
ax[2].hist2d(data[2][1],timeDiff0_4, bins = (LYSObinRange,timebinsRange),cmin = 1)
ax[3].hist2d(data[2][1],timeDiff2_4, bins = (LYSObinRange,timebinsRange),cmin = 1)

plt.show()
'''
normalize = False
LogScale = True

ax[0].hist(CsBins1[:-1],bins =CsBins1,density = normalize,weights=CsInt1/runTime,histtype = 'step',label = 'Cs-137 Ch 0', log = LogScale)
ax[1].hist(CsBins2[:-1],bins =CsBins2,density = normalize,weights=CsInt2/runTime,histtype = 'step',label = 'Cs-137 Ch 1', log = LogScale)

ax[0].hist(bckBins1[:-1],bins =bckBins1,density = normalize,weights=bckInt1/bckRunTime,histtype = 'step',label = 'Background Ch 0', log = LogScale)
ax[1].hist(bckBins2[:-1],bins =bckBins2,density = normalize,weights=bckInt2/bckRunTime,histtype = 'step',label = 'Background Ch 1', log = LogScale)
'''

    
