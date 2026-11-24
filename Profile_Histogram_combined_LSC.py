import numpy as np
import ROOT
import xml.etree.ElementTree as ET
from pathlib import Path
import scipy.stats as Stats
import pandas as pd

#This script will fit the daily calibration data with a profile histogram. 
#Note that this code does require that the file has .root format. It assumes
#that this root file has been created by CoMPASS and has been sorted into 
#coincidences from the ROOT_data_reduction.py script that I have written.


#########################################################################################
# Global Settings
#########################################################################################
# ROOT.gROOT.SetBatch(True) #Makes it so that the plots don't all show when they are created. 
# ROOT.EnableImplicitMT() #Enables implicit multithreading to speed up the code. 

ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)

ROOT.gStyle.SetLabelSize(0.05, "XYZ")   # axis tick-label size (numbers along the axis)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")   # axis title size (e.g. "Energy (ADC)")
ROOT.gStyle.SetTitleFontSize(0.06)      # the histogram's own title (top of plot)
ROOT.gStyle.SetTextSize(0.05)           # general text (legends, TLatex, etc., unless overridden)
ROOT.gStyle.SetLegendTextSize(0.1)   

#########################################################################################

class coincData:
    def __init__(self,file,nBins,LSCch,NaIch,chLabels,scaleFactor,saveFP,ROOTFile = True):
        

        self.bins = nBins
        self.scaleFactor = scaleFactor
        self.fileName = file.stem #Grab the file name
        self.chLabels = chLabels #Store the labels from the settings.xml file
        self.LSCch = LSCch #Array of LSC Channels
        self.NaIch = NaIch #Array of the single NaI channels 
        
        if "Left" in self.chLabels[f'Channel {LSCch[0]}'][0]:
            self.Lch = LSCch[0]
            self.Rch = LSCch[1]
        else:
            self.Lch = LSCch[1]
            self.Rch = LSCch[0]
        
        self.saveFP = saveFP
        
        self.cPairsStr = '_'.join([str(i) for i in self.LSCch + [self.NaIch]])
        if ROOTFile:
            self.file = ROOT.TFile.Open(str(file)) #Open the file as a root TFile
            self.cols = []
            for ch in self.LSCch + [self.NaIch]:
                ChCol = [f"Timestamp_{ch}",f"Board_{ch}",f"Energy_{ch}",f"EnergyShort_{ch}",f"Flags_{ch}"]
                self.cols += ChCol
                
            
            self.tree = self.file.Get(f"coinc_ch_{self.cPairsStr}")
            
            self.df = ROOT.RDataFrame(self.tree)
        else:
            self.file = file
            self.cols = []
            for ch in self.LSCch + [self.NaIch]:
                ChCol = [f"Timestamp_{ch}",f"Energy_{ch}"]
                self.cols += ChCol
            
            self.PdDf = pd.read_csv(self.file, header = 0, names = self.cols, usecols = [2,3,8,9,14,15], skiprows = 4, sep = ';')
            
            self.npColumns = {c: np.ascontiguousarray(self.PdDf[c].to_numpy()) for c in self.PdDf.columns}
            
            self.df = ROOT.RDF.FromNumpy(self.npColumns)
        
        self.runTime = (self.df.Max(f"Timestamp_{self.Lch}").GetValue() - self.df.Min(f"Timestamp_{self.Lch}").GetValue())*1e-12
        
    def combineLSC(self):
        colName = 'LSC_Combined'
        self.df = self.df.Define(colName, f"std::sqrt(Energy_{self.Lch} * Energy_{self.Rch})")
    
        
    def getCoincData(self):
        #Create the 1D and 2D Histograms for use in the class. 
        
        self.combineLSC()
        
        self.CombLSC1D = self.df.Histo1D(
            (f"Combined_LSC_1D_Hist_NaI_ch_{self.NaIch}",
             f"LSC Combined; Combined LSC Integral (ADC); Counts",
             self.bins,0,4000),
            "LSC_Combined"
        )
        
        self.NaI1D = self.df.Histo1D(
            (f"NaI_ch_{self.NaIch}",
             f"{self.chLabels[f'Channel {self.NaIch}'][0]}; {self.chLabels[f'Channel {self.NaIch}'][0]} Integral (ADC); Counts",
             self.bins,0,4000),
            f"Energy_{self.NaIch}"
        )
        
        self.LSCL1D = self.df.Histo1D(
            (f"LSC_left_NaI_ch_{self.NaIch}",
             f"{self.chLabels[f'Channel {self.Lch}'][0]}; {self.chLabels[f'Channel {self.Lch}'][0]} Integral (ADC); Counts",
             self.bins,0,4000),
            f"Energy_{self.Lch}"
        )
        
        self.LSCR1D = self.df.Histo1D(
            (f"LSC_right_NaI_ch_{self.NaIch}",
             f"{self.chLabels[f'Channel {self.Rch}'][0]}; {self.chLabels[f'Channel {self.Rch}'][0]} Integral (ADC); Counts",
             self.bins,0,4000),
            f"Energy_{self.Rch}"
        )
        
        self.LSCL2D = self.df.Histo2D(
            (f"LSC_Left_NaI_ch_{self.NaIch}_2D_Hist",
             f"LSC Left {self.chLabels[f'Channel {self.NaIch}'][0]} 2D Hist; {self.chLabels[f'Channel {self.NaIch}'][0]} Integral (ADC);{self.chLabels[f'Channel {self.Lch}'][0]} Integral (ADC); Counts",
             self.bins,0,4000,
             self.bins,0,4000),
            f"Energy_{self.NaIch}",
            f"Energy_{self.Lch}"
        )
        
        self.LSCR2D = self.df.Histo2D(
            (f"LSC_Right_NaI_ch_{self.NaIch}_2D_Hist",
             f"LSC Right {self.chLabels[f'Channel {self.NaIch}'][0]} 2D Hist; {self.chLabels[f'Channel {self.NaIch}'][0]} Integral (ADC);{self.chLabels[f'Channel {self.Rch}'][0]} Integral (ADC); Counts",
             self.bins,0,4000,
             self.bins,0,4000),
            f"Energy_{self.NaIch}",
            f"Energy_{self.Rch}"
        )
        
        self.CombLSC2D = self.df.Histo2D(
            (f"Combined_LSC_NaI_ch_{self.NaIch}_2D_Hist",
             f"Combined LSC {self.chLabels[f'Channel {self.NaIch}'][0]} 2D Hist; {self.chLabels[f'Channel {self.NaIch}'][0]} Integral (ADC);Combined LSC Integral (ADC); Counts",
             self.bins,0,4000,
             self.bins,0,4000),
            f"Energy_{self.NaIch}",
            "LSC_Combined"
        )
        
        self.LSC2D = self.df.Histo2D(
            (f"LSC_Right_LSC_Left_ch_{self.NaIch}_2D_Hist",
             f"LSC Right LSC Left 2D Hist; {self.chLabels[f'Channel {self.Rch}'][0]} Integral (ADC);{self.chLabels[f'Channel {self.Lch}'][0]}Integral (ADC); Counts",
             self.bins,0,4000,
             self.bins,0,4000),
            f"Energy_{self.Lch}",
            f"Energy_{self.Rch}"
        )
        
    def plotLSCData(self):
        
        self.runTime = (self.df.Max(f"Timestamp_{self.Lch}").GetValue() - self.df.Min(f"Timestamp_{self.Lch}").GetValue())*1e-12
        
        self.LSCCanvas = ROOT.TCanvas(f'LSC_Canvas_coinc_ch_{self.cPairsStr}',"LSC Plots", 1400,1400)
        self.LSCCanvas.Divide(2,2)
        
        self.LSCL1D.SetLineColor(ROOT.kBlue)
        self.LSCR1D.SetLineColor(ROOT.kRed)
        self.CombLSC1D.SetLineColor(ROOT.kBlack)
        
        for i,ch in enumerate([self.LSCL1D,self.LSCR1D,self.LSC2D,self.CombLSC1D]):
            
            pad = self.LSCCanvas.cd(i+1)
            canvasPadding(pad,0.15,0.15,0.15,0.15)
            
            h = ch.GetValue()
            h.SetStats(0)
            if h.InheritsFrom("TH2"):
                h.Draw("COLZ")
            else:
                h.Sumw2()
                h.Scale(1/self.runTime)
                h.GetYaxis().SetTitle("Counts/Bin/Second")
                h.Draw("HIST")
                
        self.LSCL1D.Draw("HIST SAME")
        self.LSCR1D.Draw("HIST SAME")    
        
        self.LSCCanvasLegend = ROOT.TLegend(0.6,0.6,0.9,0.9)
        self.LSCCanvasLegend.SetBorderSize(1)
        self.LSCCanvasLegend.SetFillColorAlpha(0,0)
        self.LSCCanvasLegend.SetTextSize(0.04)
        
        self.LSCCanvasLegend.AddEntry(self.LSCL1D.GetValue(),"LSC Left PMT",'l')
        self.LSCCanvasLegend.AddEntry(self.LSCR1D.GetValue(),"LSC Right PMT",'l')
        self.LSCCanvasLegend.AddEntry(self.CombLSC1D.GetValue(),"LSC Combined",'l')
        
        self.LSCCanvasLegend.Draw()
                
        self.LSCCanvas.SaveAs(str(self.saveFP / f'LSC_hists_channels_{self.cPairsStr}.png'))
        self.LSCCanvas.SaveAs(str(self.saveFP / f'LSC_hists_channels_{self.cPairsStr}.pdf'))
                

    def CreateProfiles(self,cutEndPoints):
        self.cut = ROOT.TCutG(f"Banna_cut_NaI_ch_{self.NaIch}",4)
        self.cut.SetLineColor(ROOT.kRed+2)
        self.cut.SetLineWidth(3)
        
        for i,p in enumerate(cutEndPoints):
            self.cut.SetPoint(i,p[0],p[1])
            
        self.cut.SetLineColor(ROOT.kRed+2)
        self.cut.SetLineWidth(3)
        
        #Define the cut histograms for the 3 histograms
        self.cutHistL = ROOT.TH2D(
            f"Cut_hist_LSC_Left_NaI_{self.NaIch}",
            f"LSC Left Cut Histogram",
            self.bins, 0,4000,self.bins,0,4000)
        
        self.cutHistR = ROOT.TH2D(
            f"Cut_hist_LSC_Right_NaI_{self.NaIch}",
            f"LSC Right Cut Histogram",
            self.bins, 0,4000,self.bins,0,4000)
        
        self.cutHistComb = ROOT.TH2D(
            f"Cut_hist_comb_LSC_NaI_{self.NaIch}",
            f"Combined LSC Cut Histogram",
            self.bins, 0,4000,self.bins,0,4000)
        
        #Define the profiles for the three histograms. 
        self.profileL = ROOT.TProfile(
            f"Profile_hist_LSC_Left_NaI_ch_{self.NaIch}", "Profile Histogram",
            self.bins,0,4000)
        
        self.profileR = ROOT.TProfile(
            f"Profile_hist_LSC_Right_NaI_ch_{self.NaIch}", "Profile Histogram",
            self.bins,0,4000)
                
        self.profileComb = ROOT.TProfile(
            f"Profile_hist_Combined_LSC_NaI_ch_{self.NaIch}", "Profile Histogram",
            self.bins,0,4000)
        
        maxBinL = findMaxBin(self.LSCL2D,self.cut)
        maxBinR = findMaxBin(self.LSCR2D,self.cut)
        maxBinComb = findMaxBin(self.CombLSC2D,self.cut)
        
        
        self.cutHistL,self.profileL = fillCutHist(self.LSCL2D,self.cutHistL,self.profileL,maxBinL*self.scaleFactor,self.cut)
        
        self.cutHistR,self.profileR = fillCutHist(self.LSCR2D,self.cutHistR,self.profileR,maxBinR*self.scaleFactor,self.cut)
        
        self.cutHistComb,self.profileComb = fillCutHist(self.CombLSC2D,self.cutHistComb,self.profileComb,maxBinComb*(self.scaleFactor-0.2),self.cut)

        
        self.xRangeL = findProfileRange(self.profileL)
        self.xRangeR = findProfileRange(self.profileR)
        self.xRangeComb = findProfileRange(self.profileComb)
        
        self.fitL = ROOT.TF1('Lin_Fit_L', "pol1", self.xRangeL[0],self.xRangeL[1])
        self.fitR = ROOT.TF1('Lin_Fit_R', "pol1", self.xRangeR[0],self.xRangeR[1])
        self.fitComb = ROOT.TF1('Lin_Fit_Comb', "pol1", self.xRangeComb[0],self.xRangeComb[1])
        
        for f in [self.fitL,self.fitR,self.fitComb]:
            f.SetLineColor(ROOT.kRed)
            f.SetLineWidth(3)
            
        
        self.profileL.Fit(self.fitL,"RN")
        self.profileR.Fit(self.fitR,"RN")
        self.profileComb.Fit(self.fitComb,"RN")
    
        
    def plotProfiles(self):
        self.profileCanvas = ROOT.TCanvas(f"Profile_Canvas_NaI_{self.NaIch}", "Profile Canvas", 2800,800)
        
        self.profileCanvas.Divide(3,1)
        
        profiles = [self.profileL,self.profileR,self.profileComb]
        Hists2D = [self.LSCL2D,self.LSCR2D,self.CombLSC2D]
        fits = [self.fitL,self.fitR,self.fitComb]
        labels = ["LSC Left","LSC Right", "Combined LSC"]
        self.legends = []
        self.p0,self.p1,self.p0Err,self.p1Err = [],[],[],[]
        self.chi2,self.ndof,self.pVal = [],[],[]
        
        for i, data in enumerate(zip(profiles,Hists2D,fits)):
        
            pad = self.profileCanvas.cd(i+1)
            canvasPadding(pad,0.2,0.25,0.1,0.35)
            
            data[1].SetStats(0)
            data[1].GetXaxis().SetNdivisions(505)
            data[1].GetZaxis().SetTitleOffset(1.4)
            data[1].Draw("COLZ")
            
            leg = ROOT.TLegend(0.03,0.02,0.97,0.25)
            # leg.SetNColumns(2)
            leg.SetBorderSize(1)
            leg.SetFillColorAlpha(0,0)
            leg.SetTextSize(0.05)

            self.p0.append(data[2].GetParameter(0))
            self.p1.append(data[2].GetParameter(1))
            
            self.p0Err.append(data[2].GetParError(0))
            self.p1Err.append(data[2].GetParError(1))
            
            self.chi2.append(data[2].GetChisquare())
            self.ndof.append(data[2].GetNDF())
            
            self.pVal.append(Stats.chi2.sf(self.chi2[-1], self.ndof[-1]))
            
            leg.AddEntry(data[2],f"p0 = {self.p0[-1]:.2f} #pm {self.p0Err[-1]:.2f}",'l')
            leg.AddEntry((0),f"p1 = {self.p1[-1]:.2f} #pm {self.p1Err[-1]:.2f}",'')
            
            leg.AddEntry((0),f"#chi^{{2}} / ndof = {self.chi2[-1]:.2f} / {self.ndof[-1]}",'')
            leg.AddEntry(data[2],f"p-Value = {self.pVal[-1]:.4f}",'')
            
            self.legends.append(leg)
            self.legends[-1].Draw()
            
            
            
            data[0].SetStats(0)
            data[0].SetLineColor(ROOT.kBlack)
            data[0].SetLineWidth(3)
            data[0].Draw("SAME")
            
            
            
            data[2].Draw("Same")
            
            self.cut.Draw("SAME")
            
            self.profileCanvas.Update()
            
            
        self.profileCanvas.SaveAs(str(self.saveFP / f'Profile_histogram_channels_{self.cPairsStr}.png'))
        self.profileCanvas.SaveAs(str(self.saveFP / f'Profile_histogram_channels_{self.cPairsStr}.pdf'))
            
        
        
        
        
        

        
def findProfileRange(profile,minValue = 1):
    
    ax = profile.GetXaxis()
    xBins = [x for x in range(1,profile.GetNbinsX()+1) if profile.GetBinEntries(x) >= minValue]
    
    return [ax.GetBinLowEdge(xBins[0]), ax.GetBinUpEdge(xBins[-1])]
    
        
        
def findMaxBin(hist,cut):
    #Small function that loops through the data and finds what bin has the most counts. 
    
    maxBin = 0
    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsY()+1):
            x = hist.GetXaxis().GetBinCenter(binx)
            y = hist.GetYaxis().GetBinCenter(biny)
            z = hist.GetBinContent(binx,biny)
            
            if not cut.IsInside(x,y) and z >0:
                if z > maxBin:
                    maxBin = z
                    
    return maxBin

def fillCutHist(hist,cuthist,profile,threshold,cut):
    
    for binx in range(1,hist.GetNbinsX()+1):
        for biny in range(1,hist.GetNbinsY()+1):
            x = hist.GetXaxis().GetBinCenter(binx)
            y = hist.GetYaxis().GetBinCenter(biny)
            z = hist.GetBinContent(binx,biny)
            
            if not cut.IsInside(x,y) and z >threshold:
                profile.Fill(x,y,z)
                cuthist.Fill(x,y,z)
                    
    return cuthist,profile
        
def canvasPadding(pad,L,R,T,B):
    pad.SetLeftMargin(L)
    pad.SetRightMargin(R)   # extra room for COLZ palette
    pad.SetTopMargin(T)
    pad.SetBottomMargin(B)       




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


#File path for the file to be fit. 
# filePath = Path("/home/nick/PhD/KDK+/Pilot_exp/2026_09_18/2026_09_18_Pilot_experiment_cal_data/RAW/coinc_sorted/SDataR_2026_09_18_Pilot_experiment_cal_data_sorted_coinc.root")

filePath = Path("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_06_19/2026_06_19_Daily_LSC_calibration_Cs137_coinc/RAW/coinc_sorted_500ns/")




#Define the output directory and create it if it doesn't exist. 


################################################################################
scaleFactor = 0.6 #Sets the scale factor that creates the cut histogram based on number of counts. 

#Currently I have you adjust the settings here to define what channels correspond to the NaI and the LSC. This potentially could be updated later based on the channel labels. 

# LSCChannels = [0,1] #Define the channel numbers for the LSC
# NaIChannels = [2,3,4,5] #Define the channel numbers for the NaI

LSCChannels = [4,5] #Define the channel numbers for the LSC
NaIChannels = [8,10,12,14] #Define the channel numbers for the NaI

ROOTFile = False

nbins = 100
combineLSC = True #Set this to true for combine LSC vs NaI. 
#################################################################################

if ROOTFile:
    #Grabs the settings.xml file for the data file provided. 
    settingsfilepath = filePath.parent.parent.parent / "settings.xml"  
    saveFilePath = filePath.parent.parent / 'Daily_calibration_fits' / 'figures'

else:
    #Grabs the settings.xml file for the data file provided. 
    settingsfilepath = filePath.parent.parent / "settings.xml"
    saveFilePath = filePath.parent / 'Daily_calibration_fits' / 'figures'
    

saveFilePath.mkdir(parents = True, exist_ok = True)
#Read in the settings.xml file to grab the channel labels defined in CoMPASS. 
chNames = ReadInChannelNames(settingsfilepath)


#Define the cut endpoints for a rough cut around the shape of the histogram. 

# cutEndPoints = [[0,0], #[x coords, y coords]
#                 [0,4000],
#                 [50,4000],
#                 [50,2000],
#                 [1200,50],
#                 [4000,50],
#                 [4000,0],
#                 [0,0]]


cutEndPoints = [[0,0], #[x coords, y coords]
                [0,4000],
                [50,4000],
                [50,800],
                [800,50],
                [4000,50],
                [4000,0],
                [0,0]]


data = []
for NaIch in NaIChannels:
    if ROOTFile:
        cData = coincData(filePath,nbins, LSCChannels,NaIch,chNames,scaleFactor,saveFilePath,ROOTFile)
    else:
        fileName = filePath.parent.parent.stem
        print(f"FileName: {fileName}")
        # file = filePath / f"SDataR_2026_04_28_Daily_LSC_calibration_Cs137_coinc_2_coinc_{LSCChannels[0]}_{LSCChannels[1]}_{NaIch}.txt"
        
        file = filePath / f"SDataR_{fileName}_coinc_{LSCChannels[0]}_{LSCChannels[1]}_{NaIch}.txt"
        
        cData = coincData(file,nbins, LSCChannels,NaIch,chNames,scaleFactor,saveFilePath,ROOTFile)
        
    data.append(cData)
    cData.getCoincData()
    cData.plotLSCData()
    cData.CreateProfiles(cutEndPoints)
    cData.plotProfiles()
    
    
#Save the fit results to a txt file

saveFile = saveFilePath.parent / "Profile_Fit_Results.txt"
date = filePath.parent.parent.parent.parent.stem

with open(saveFile, "w") as f:
    f.write("Date,NaI_ch,")
    f.write("Left_po,Left_po_err,Left_p1,Left_p1_err,Left_chi2,left_ndof,left_pval,")
    f.write("Right_po,Right_po_err,Right_p1,Right_p1_err,Right_chi2,Right_ndof,Right_pval,")
    f.write("Comb_po,Comb_po_err,Comb_p1,Comb_p1_err,Comb_chi2,Comb_ndof,Comb_pval\n")

    for d in data:
        f.write(f"{date},{d.NaIch},")
        for i in range(3):
            f.write(f"{d.p0[i]},{d.p0Err[i]},{d.p1[i]},{d.p1Err[i]},{d.chi2[i]},{d.ndof[i]},{d.pVal[i]}")
            if i != 2:
                f.write(',')
        f.write("\n")
        
    f.close()


