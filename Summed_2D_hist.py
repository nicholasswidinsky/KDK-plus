import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import time 
import datetime


startTime = time.time()

class ChData:
    def __init__(self,ch,E):
        self.ch = int(ch)
        self.E = [E]
        # self.t = [t]
        
    def AddEvent(self,E):
        self.E.append(E)
        # self.t.append(t)
    
        
        
class coincData:
    def __init__(self,channel,Energy):
        self.chList = channel
        self.chData = {}
        for i,ch in enumerate(channel):
            self.chData.update({f'{ch}' : ChData(int(ch),Energy[i])})
        
        # print(self.chData['4'].ch)
        
    def AddEvent(self,channel,Energy):
        
        for i, ch in enumerate(channel):
            self.chData[f'{ch}'].AddEvent(Energy[i])
    

class Summed2DHist:
    def __init__(self,ch):
        self.dates = []
        self.Channels = ch #List of the channels that are being summed (LSC is always along the x-axis. NaI are always along the y-axis)
        self.xInts = np.array([]) #Integrals along the x axis
        self.yInts = np.array([]) #Integrals along the y axis
        
        self.SeparationPoints = []
        
    def addHistData(self,xData,yData, date):
        self.xInts = np.concatenate((self.xInts,xData)) #Adds the latest data to the x and y arrays
        self.yInts = np.concatenate((self.yInts,yData))
        
        if date in self.dates:
            self.SeparationPoints[-1] + len(self.xInts)
        else:
            self.dates.append(date)
            self.SeparationPoints.append(len(self.xInts)) #Adds the separation point to note where each day stops. 
        
    
    def plot2DHist(self,savefilepath):
        
        bins = 100
        binRange = [0,4000]
        fig,ax = plt.subplots(6,6,figsize = (40,40))
        plt.tight_layout()
        plt.subplots_adjust(left = 0.05, bottom = 0.05,hspace = 0.2,wspace = 0.2)
        axsFlat = ax.flatten()
        
        for i,axs in enumerate(axsFlat):
            axs.hist2d(self.xInts[0:self.SeparationPoints[i]],self.yInts[0:self.SeparationPoints[i]], bins = [bins,bins], range = [binRange,binRange], cmin = 1)
            axs.set_xlabel('LSC Int (ADC)')
            axs.set_ylabel('NaI Int (ADC)')
            axs.set_title(f'Day {i+1}')
            
        savepath = savefilepath / 'Summed_Histograms'
        savepath.mkdir(parents=True,exist_ok=True)
        plt.savefig(savepath / f'Summed_hist_ch_{self.Channels[0]}_{self.Channels[1]}.png')  
        plt.savefig(savepath / f'Summed_hist_ch_{self.Channels[0]}_{self.Channels[1]}.pdf')  
        

def readInData(file):
    with open(file) as f: 
        for i in range(3): #skips the first 3 lines of the file. 
            next(f)
            
        for i, line in enumerate(f):
            if i == 0: #Look at the third line to get the total number of channels in this coincidence. 
                numChannels = int(line.split("\n")[0].split(": ")[1])
            else:
                data = line.split("\n")[0].split(";")
                ch, E, t = [],[],[]
                for j in range(numChannels):
                    ch.append(int(data[6*j+1]))
                    # t.append(float(data[6*j+2])/1e3)
                    E.append(float(data[6*j+3]))
                if i == 1: #Use the first line to initialize the class objects using the channels. 
                    cData = coincData(ch,E)
                else:
                    cData.AddEvent(ch,E)
                    
    return cData

# def readInData(file):
#     with open(file) as f:
#         lines = f.readlines()

#     numChannels = int(lines[3].split(": ")[1])

#     # Parse all data rows at once
#     data_lines = [line.strip().split(";") for line in lines[4:] if line.strip()]
#     raw = np.array(data_lines)

#     # Extract all channels, times, energies as vectors
#     ch_cols = [raw[:, 6*j + 1] for j in range(numChannels)]
#     t_cols  = [raw[:, 6*j + 2].astype(float) / 1e3 for j in range(numChannels)]
#     E_cols  = [raw[:, 6*j + 3].astype(float) for j in range(numChannels)]

#     channels = [int(col[0]) for col in ch_cols]

#     # Initialize with first event (preserves your existing constructor signature)
#     cData = coincData(channels, [t[0] for t in t_cols], [E[0] for E in E_cols])

#     # Bulk-load remaining events into ChData, bypassing the slow append loop
#     for ch, E_vec, t_vec in zip(channels, E_cols, t_cols):
#         cData.chData[ch].E = list(E_vec)
#         cData.chData[ch].t = list(t_vec)

#     return cData

    
    
def ReadInChannelNames(settings):
    ##############################################################################################
    #   Function that reads in the settings.xml file produced by CoMPASS to grab the channel     #
    #   labels and make a dictionary for them. Note that I did use chatGPT to make this function #
    #   As much as it shames me, I had no idea how to read in this file so I cheated a bit.      #
    #   I have gone through the code, and the comments that ChatGPT left and it all makes sense. #
    #   I have also tested it to ensure that it is working properly.                             #
    ##############################################################################################    

    # Parse the XML file into a tree structure
    tree = ET.parse(settings)
    root = tree.getroot()

    # Initialize an empty dictionary to store channel index and label pairs
    channels = {}

    # Iterate over every <channel> element in the file
    for ch in root.iter("channel"):
        # Find the <index> tag within this <channel> element
        index_elem = ch.find("index")
        
        # Find the <values> tag, which contains parameter entries like labels
        values_elem = ch.find("values")
        
        # Continue only if both index and values are found
        if index_elem is not None and values_elem is not None:
            label = None  # Placeholder for the channel label text

            # Search through each <entry> element under <values>
            for entry in values_elem.findall("entry"):
                # Each entry has a <key> and <value>
                key_elem = entry.find("key")
                value_elem = entry.find("value")
                
                # Check if this entry is the label for the channel
                if key_elem is not None and key_elem.text == "SW_PARAMETER_CH_LABEL":
                    # If found, extract the label text (safely)
                    label = value_elem.text.strip() if value_elem is not None else None
                    break  # No need to check further entries for this channel

            # If a label was found, store it in the dictionary with its index
            if label is not None:
                channels[f'Channel {int(index_elem.text)}'] = [label]

    # Print the resulting dictionary: {channel_index: label, ...}
    return channels


rootfilePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/') #filepath to the coinc sorted directory. 

LSCChannels = [4,5]
NaIChannels = [8,10,12,14] #Protects against the possibility of having a werid channel coincidence layout with incorrect channel numbers. 

# LSCChannels = [0,1]
# NaIChannels = [2,3,4,5]


pattern = re.compile(rf'_coinc_{LSCChannels[0]}_{LSCChannels[1]}_({NaIChannels[0]}|{NaIChannels[1]}|{NaIChannels[2]}|{NaIChannels[3]}).txt$')

coincFiles = sorted([
    f for f in rootfilePath.glob(f'**/*_coinc_{LSCChannels[0]}_{LSCChannels[1]}_*.txt')
    if pattern.search(f.name)
])
histogramData = []
for i in LSCChannels:
    for j in NaIChannels:
        histogramData.append(Summed2DHist([i,j]))

readInStartTime = time.time()
print('Reading in histogram data:')

for i,filePath in enumerate(coincFiles):
    date = filePath.parent.parent.parent.parent.stem
    settingsFilePath = filePath.parent.parent.parent / 'settings.xml'
    # saveFilePath = filePath.parent.parent / 'Daily_calibration_fits' / 'figures'
    # outputFileName = filePath.stem
    # outputFile = saveFilePath / f"{outputFileName}_Profile_hist_fit_results.txt"
    
    # if outputFile.exists():
    #     outputFile.unlink()

    # saveFilePath.mkdir(parents=True, exist_ok=True)
    savefilepath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/')
    if i == 0:
        Detectors = ReadInChannelNames(settingsFilePath)
    cData = readInData(filePath)
    chPairs = [[cData.chList[0],cData.chList[2]], [cData.chList[1],cData.chList[2]]] #Makes a 2D array with the two channel pairs for the summed histogram data. 
    print(cData.chList)


    for j,hist in enumerate(histogramData):
        if hist.Channels == chPairs[0]:
            hist.addHistData(cData.chData[str(chPairs[0][0])].E,cData.chData[str(chPairs[0][1])].E,date)
        elif hist.Channels == chPairs[1]:
            hist.addHistData(cData.chData[str(chPairs[1][0])].E,cData.chData[str(chPairs[1][1])].E,date)


totalReadTime = time.time()
print(f'Total time to read in data: \t {totalReadTime - readInStartTime} s')

for hist in histogramData:
    print(f'Plotting Histogram from channels: {hist.Channels}')
    hist.plot2DHist(savefilepath)



    
    
    
    