import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from pathlib import Path
import scipy.stats as Stats
from ROOT import TCanvas, TH2D, TCutG,TProfile, TF1, kRed,TLegend


class ChData:
    def __init__(self,ch,E,t):
        self.ch = ch
        self.E = [E]
        self.t = [t]
        
    def AddEvent(self,E,t):
        self.E.append(E)
        self.t.append(t)
        
        
class coincData:
    def __init__(self,channel,time,Energy):
        self.chList = channel
        self.chData = {}
        for i,ch in enumerate(channel):
            self.chData.update({ch : ChData(int(ch),Energy[i],time[i])})
        
        # print(self.chData['4'].ch)
        
    def AddEvent(self,channel,time,Energy):
        
        for i, ch in enumerate(channel):
            self.chData[ch].AddEvent(Energy[i],time[i])
            
            
class Hists2D:
    def __init__(self,ch,Energy):
        self.ch = ch
        self.E = Energy
        self.data = {f"{ch[0]}":Energy[0],
                     f"{ch[1]}":Energy[1]}
        self.hist2D = TH2D(f"{channels[f'Channel {self.ch[0]}'][0]} vs {channels[f'Channel {self.ch[1]}'][0]} 2D Histogram", f"{channels[f'Channel {self.ch[0]}'][0]} vs {channels[f'Channel {self.ch[1]}'][0]}", 100, 0, 4000, 100, 0, 4000)
        for i in range(len(Energy[0])):
            self.hist2D.Fill(Energy[0][i],Energy[1][i])
        self.profile = None
        self.cut = None  
        self.fit = None      
        
    def Draw2DHist(self):
        self.canvas = TCanvas()
        
        self.hist2D.GetXaxis().SetTitle(f"{channels[f'Channel {self.ch[0]}'][0]} Energy (ADU)")
        self.hist2D.GetYaxis().SetTitle(f"{channels[f'Channel {self.ch[1]}'][0]} Energy (ADU)")
        self.hist2D.GetZaxis().SetTitle("Counts")
        
        self.hist2D.SetStats(0)
        self.hist2D.Draw("COLZ")
        if self.cut is not None:
            self.cut.Draw("same")
        if self.profile is not None:
            self.profile.SetStats(0) 
            self.profile.Draw("same")
        if self.fit is not None:
            self.fit.Draw("same")
            
            # Extract fit results
            p0       = self.fit.GetParameter(0)
            p0_err   = self.fit.GetParError(0)
            p1       = self.fit.GetParameter(1)
            p1_err   = self.fit.GetParError(1)
            chi2     = self.fit.GetChisquare()
            ndof     = self.fit.GetNDF()
            p_value  = Stats.chi2.sf(chi2, ndof)  # survival function = 1 - CDF
            
            # Build legend — positioned to the left of the colour bar
            # (x1, y1, x2, y2) in NDC coordinates (0-1)
            legend = TLegend(0.12, 0.65, 0.55, 0.88)
            legend.SetHeader("Fit Results (pol1)", "C")
            legend.AddEntry(self.fit, f"p0 = {p0:.2f} #pm {p0_err:.2f}", "l")
            legend.AddEntry((0), f"p1 = {p1:.4f} #pm {p1_err:.4f}", "")
            legend.AddEntry((0), f"#chi^{{2}} / ndof = {chi2:.2f} / {ndof}", "")
            legend.AddEntry((0), f"p-value = {p_value:.4f}", "")
            legend.SetBorderSize(1)
            legend.SetFillColorAlpha(0, 0.6)  # semi-transparent background
            self.legend = legend               # store to prevent garbage collection
            self.legend.Draw()
            
            
        
        self.canvas.Update()
        
    def cutHist(self,cutEndPoints):
        self.cut = TCutG(f'Banana_cut_{channels[f'Channel {self.ch[0]}'][0]}_vs_{channels[f'Channel {self.ch[1]}'][0]}', 4)
        for i in range(len(cutEndPoints)):
            self.cut.SetPoint(i,cutEndPoints[i][0], cutEndPoints[i][1])
        self.cut.SetLineColor(kRed+2)
        self.cut.SetLineWidth(3)
        
    def makeProfileHist(self,fitRange):
        self.profile = TProfile(f"profile_ch_{channels[f'Channel {self.ch[0]}'][0]}_vs_ch_{channels[f'Channel {self.ch[1]}'][0]}",f"profile Histogram {channels[f'Channel {self.ch[0]}'][0]} vs {channels[f'Channel {self.ch[1]}'][0]}", 100,0,4000,0,4000,"")
        for bini in range(1, self.hist2D.GetNbinsX() + 1):
            for binj in range(1, self.hist2D.GetNbinsY() + 1): 
                x = self.hist2D.GetXaxis().GetBinCenter(bini)
                y = self.hist2D.GetYaxis().GetBinCenter(binj)
                z = self.hist2D.GetBinContent(bini,binj)
                
                if not self.cut.IsInside(x,y) and z > 0:
                    self.profile.Fill(x,y,z)
        self.fit = TF1(f"fit_{channels[f'Channel {self.ch[0]}'][0]}_vs_{channels[f'Channel {self.ch[1]}'][0]}","pol1",fitRange[0], fitRange[1])
        self.profile.Fit(self.fit, "RN")
        self.fit.SetLineColor(kRed)
                
                
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
                    ch.append(data[6*j+1])
                    t.append(float(data[6*j+2])/1e3)
                    E.append(float(data[6*j+3]))
                if i == 1: #Use the first line to initialize the class objects using the channels. 
                    cData = coincData(ch,t,E)
                else:
                    cData.AddEvent(ch,t,E)
                    
    return cData
  
  
def make2DHists(cData):
    
    hist2DData = []
    for i,ch in enumerate(cData.chList):
        j = i+1
        while j < len(cData.chList):
            hist2DData.append(Hists2D([cData.chList[i], cData.chList[j]],
                                      [cData.chData[f'{cData.chList[i]}'].E, cData.chData[f'{cData.chList[j]}'].E]))
            j+=1
    
    return hist2DData

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
    
    
# filePath = '/home/nick/PhD/KDK+/Annulus_Compton_scatter_V1/2026_01_28/NaI_annulus_LS_2_Cs137_All_NaI_higher_HV/RAW/coinc_sorted/SDataR_NaI_annulus_LS_2_Cs137_All_NaI_higher_HV_coinc_4_5_8.txt'

filePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_04_30/2026_04_30_Daily_LSC_calibration_Cs137_coinc/RAW/coinc_sorted/SDataR_2026_04_30_Daily_LSC_calibration_Cs137_coinc_coinc_4_5_8.txt')

settingsFilePath = filePath.parent.parent.parent / 'settings.xml'

#Read in the detector names from the settings file. 
channels = ReadInChannelNames(settingsFilePath)
                
cData = readInData(filePath)

hist2DData = make2DHists(cData)

cutEndPoints = [[0,0],
                [0,1000],
                [1600,0],
                [0,0]]

fitRange = [200,1100]

for data in hist2DData:
    data.cutHist(cutEndPoints)
    data.makeProfileHist(fitRange)
    data.Draw2DHist()

# test2DHist = Hists2D([cData.chData['4'],cData.chData['8']],[cData.chData['4'].E,cData.chData['8'].E])
# test2DHist.Draw2DHist()





# test2Dhist = TH2D("test","test Hist", 100,0,4000,100,0,4000)
# for i in range(len(cData.chData['4'].E)):
#     test2Dhist.Fill(cData.chData['4'].E[i],cData.chData['5'].E[i])
# test2Dhist.Draw()