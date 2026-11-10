import numpy as np
import ROOT
from pathlib import Path
import xml.etree.ElementTree as ET
from itertools import combinations
import math


#########################################################################################
# Global Settings
#########################################################################################
# ROOT.gROOT.SetBatch(True)
# ROOT.EnableImplicitMT() #Enables implicit multithreading to speed up the code. 

ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)

ROOT.gStyle.SetLabelSize(0.05, "XYZ")   # axis tick-label size (numbers along the axis)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")   # axis title size (e.g. "Energy (ADC)")
ROOT.gStyle.SetTitleFontSize(0.06)      # the histogram's own title (top of plot)
ROOT.gStyle.SetTextSize(0.05)           # general text (legends, TLatex, etc., unless overridden)
ROOT.gStyle.SetLegendTextSize(0.1)   
#########################################################################################

class positionData:
    def __init__(self,file,channels,ncols,fitWindow,saveFP):
        self.file = ROOT.TFile.Open(str(file))
        # self.file.ls() #Prints the contents of the file 
        
        settingsFilePath = file.parent.parent.parent / 'settings.xml'
        
        self.saveFP = saveFP / "Position_Fits"
        self.saveFP.mkdir(parents = True, exist_ok = True)

        
        self.chLabels = ReadInChannelNames(settingsFilePath)
        
        self.filename = file.stem
        self.positionSTR = self.filename.split("_")[6]
        
        self.channels = channels
        
        if 'R' in self.positionSTR:
            self.position = int(self.positionSTR[:-1])
        elif 'L' in self.positionSTR:
            self.position = -1 * int(self.positionSTR[:-1])
        else:
            self.position = int(self.positionSTR)

        if "Left" in self.chLabels[f'Channel 0'][0]:
            self.Lch = 0
            self.Rch = 1
        else:
            self.Lch = 1
            self.Rch = 0
   
        self.NaIch = 2
        
        
        self.fitWindow = fitWindow
        
        ################################################################
        # Define the fits 
        ################################################################
        
        fitNames = ['Ngauss', 'mu', 'sigma','Nexp', 'tau', 'erfStep']
        
        #Left PMT fits
        #No erf
        self.fitFuncL = ROOT.TF1(
            f"Position_{self.position}_fit_func_L", 'gaus(0) + expo(3)',
            fitWindow[f'{self.position}'][0][0],
            fitWindow[f'{self.position}'][0][1])
        
        
        
        self.fitFuncL.SetLineColor(ROOT.kRed)
        self.fitFuncL.SetLineWidth(3)
        
        #erf
        self.fitFuncLerf = ROOT.TF1(
            f"Position_{self.position}_fit_func_L", 'gaus(0) + expo(3) + [5]*TMath::Erfc((x-[1])/[2])',
            fitWindow[f'{self.position}'][0][0],
            fitWindow[f'{self.position}'][0][1])
        
        self.fitFuncLerf.SetLineColor(ROOT.kMagenta)
        self.fitFuncLerf.SetLineWidth(3)
        
        #Right PMT fits
        #No erf
        self.fitFuncR = ROOT.TF1(
            f"Position_{self.position}_fit_func_R", 'gaus(0) + expo(3)',
            fitWindow[f'{self.position}'][1][0],
            fitWindow[f'{self.position}'][1][1])
        
        self.fitFuncR.SetLineColor(ROOT.kRed)
        self.fitFuncR.SetLineWidth(3)
        
        #erf
        self.fitFuncRerf = ROOT.TF1(
            f"Position_{self.position}_fit_func_R", 'gaus(0) + expo(3) + [5]*TMath::Erfc((x-[1])/[2])',
            fitWindow[f'{self.position}'][1][0],
            fitWindow[f'{self.position}'][1][1])
        
        self.fitFuncRerf.SetLineColor(ROOT.kMagenta)
        self.fitFuncRerf.SetLineWidth(3)
        
        
        for i,n in enumerate(fitNames):
            self.fitFuncL.SetParName(i,n)
            self.fitFuncLerf.SetParName(i,n)
            self.fitFuncR.SetParName(i,n)
            self.fitFuncRerf.SetParName(i,n)
        
        
        self.LpadInd = self.position + ncols//2 +1 #Sets the position on the pad for the left PMT
        self.RpadInd = self.position + ncols //2 + ncols +1 #Sets the position on the pad for the right PMT
        
        self.cols = []
        for ch in channels:
            ChCol = [f"Timestamp_{ch}",f"Board_{ch}",f"Energy_{ch}",f"EnergyShort_{ch}",f"Flags_{ch}"]
            self.cols += ChCol
            
        self.cPairsStr = '_'.join([str(i) for i in channels])
        self.tree = self.file.Get(f"coinc_ch_{self.cPairsStr}")
        
        self.df = ROOT.RDataFrame(self.tree)
        
    def plot1DHist(self,canvas):
        
        #####################################################################
        # Create the two histograms 
        #####################################################################
        
        print(f"\n Plotting Position {self.positionSTR} \n")
        
        
        #Left PMT Hist
        self.LSCLHist1D = self.df.Histo1D((f"Energy_Hist_{self.Lch}", f"{self.chLabels[f'Channel {self.Lch}'][0]} Position {self.positionSTR};{self.chLabels[f'Channel {self.Lch}'][0]} Integral (ADC);Counts", 100,0,4000),f"Energy_{self.Lch}")
        
        #Right PMT Hist
        self.LSCRHist1D = self.df.Histo1D((f"Energy_Hist_{self.Rch}", f"{self.chLabels[f'Channel {self.Rch}'][0]} Position {self.positionSTR};{self.chLabels[f'Channel {self.Rch}'][0]} Integral (ADC);Counts", 100,0,4000),f"Energy_{self.Rch}")
        
        #Turn off stats
        self.LSCLHist1D.SetStats(0)
        self.LSCRHist1D.SetStats(0)
        
        
        #Swap to the first pad to plot the Left PMT data
        pad = canvas.cd(self.LpadInd)
        
        canvasPadding(pad)
        self.LSCLHist1D.Draw()
        
        #Set initial parameters for the fits. 
        ampGuess = 500
        minRange, maxRange = self.fitWindow[f'{self.position}'][0]
        sigGuess = (maxRange-minRange) / 10
        muGuess = (maxRange+minRange)/2
        
        expAmpGuess = 5
        tau = 1e-10
        stepHeight = 30
        self.fitFuncL.SetParameters(ampGuess,muGuess,sigGuess,expAmpGuess,tau)
        self.fitFuncR.SetParameters(ampGuess,muGuess,sigGuess,expAmpGuess,tau)
        
        for fitFunc in [self.fitFuncL,self.fitFuncR]:
            fitFunc.SetParLimits(1, minRange, maxRange) # mean stays in the window
            fitFunc.SetParLimits(2, 1, (maxRange - minRange)) # sigma positive
            # fitFunc.SetParLimits(4, 0, 10) # tau positive

        
        self.fitFuncLerf.SetParameters(ampGuess,muGuess,sigGuess,expAmpGuess,tau,stepHeight)
        self.fitFuncRerf.SetParameters(ampGuess,muGuess,sigGuess,expAmpGuess,tau,stepHeight)
        
        for fitFunc in [self.fitFuncLerf,self.fitFuncRerf]:
            fitFunc.SetParLimits(1, minRange, maxRange)          # mean stays in the window
            fitFunc.SetParLimits(2, 1, (maxRange - minRange))    # sigma positive
            fitFunc.SetParLimits(5, 0, 1e4)          # step height positive
            fitFunc.SetParLimits(4, 0, 10) # tau positive
            

        
        #Fit the data no erf
        self.LSCLHist1D.Fit(self.fitFuncL, "RS")
        self.fitFuncL.Draw("SAME")
        
        #Fit the data erf
        self.LSCLHist1D.Fit(self.fitFuncLerf, "RS")
        self.fitFuncLerf.Draw("SAME")
        self.LSCLHist1D.GetXaxis().SetNdivisions(505) 
       
        
        canvas.Update()
        
        pad = canvas.cd(self.RpadInd)
        canvasPadding(pad)
        self.LSCRHist1D.Draw()
        
        #Fit the data no erf
        self.LSCRHist1D.Fit(self.fitFuncR, "RS")
        self.fitFuncR.Draw("SAME")
        
        #Fit the data erf
        self.LSCRHist1D.Fit(self.fitFuncRerf, "RS")
        self.fitFuncRerf.Draw("SAME")
        self.LSCRHist1D.GetXaxis().SetNdivisions(505) 
        
        canvas.Update()
        
        
        #Create the plot for the individual position fits
        self.FitCanvas = ROOT.TCanvas(f"Individual_Fit_canvas_{self.positionSTR}", f"Position {self.positionSTR} fit", 800,1600)
        
        self.FitCanvas.Divide(1,2)
        
        pad = self.FitCanvas.cd(1)
        canvasPadding(pad)
        
        leftMin = self.fitFuncL.GetXmin()
        leftMax = self.fitFuncL.GetXmax()
        
        self.LSCLHist1D.Draw()
        self.fitFuncL.Draw("SAME")
        self.fitFuncLerf.Draw("SAME")
        
        ROOT.gPad.Update()
        Leftymin = ROOT.gPad.GetUymin()
        Leftymax = ROOT.gPad.GetUymax()
        
        self.LeftMinLine = ROOT.TLine(leftMin,Leftymin,leftMin,Leftymax)
        self.LeftMaxLine = ROOT.TLine(leftMax,Leftymin,leftMax,Leftymax)
        
        for line in [self.LeftMinLine, self.LeftMaxLine]:
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("same")
        
        self.FitCanvas.Update()
        
        pad = self.FitCanvas.cd(2)
        canvasPadding(pad)
        
       
        
        self.LSCRHist1D.Draw()
        self.fitFuncR.Draw("SAME")
        self.fitFuncRerf.Draw("SAME")
        
        rightMin = self.fitFuncR.GetXmin()
        rightMax = self.fitFuncR.GetXmax()
        
        ROOT.gPad.Update()
        
        Rightymin = ROOT.gPad.GetUymin()
        Rightymax = ROOT.gPad.GetUymax()
        
        self.RightMinLine = ROOT.TLine(rightMin,Rightymin,rightMin,Rightymax)
        self.RightMaxLine = ROOT.TLine(rightMax,Rightymin,rightMax,Rightymax)
        
        for line in [self.RightMinLine, self.RightMaxLine]:
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("same")
        
        self.FitCanvas.Update()
        
        self.FitCanvas.SaveAs(str(self.saveFP / f"Individual_fits_{self.positionSTR}.pdf"))
        self.FitCanvas.SaveAs(str(self.saveFP / f"Individual_fits_{self.positionSTR}.png"))
        

# def Fit1DData(df):

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

def ReadInNaIScale(file):
    with open(file) as f:
        next(f)
        
        for line in f:
            NaIScale = float(line.split(',')[4])
            
    return NaIScale 

def scaleNaIIntegrals(df,scale):
    df = df.Define("Energy_2_keV", f"Energy_2 / {scale}")
    
    return df

def canvasPadding(pad):
    pad.SetLeftMargin(0.2)
    pad.SetRightMargin(0.10)   # extra room for COLZ palette
    pad.SetTopMargin(0.10)
    pad.SetBottomMargin(0.15)
    

RootDirectory = Path("/home/nick/PhD/KDK+/Large_LSC_testing/Position_test_no_outer_sleeve/")
NaIScaleFile = Path("/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/Loaded_BNL_L_LSC_Small_NaI_LY_coinc_2/RAW/SDataR_Loaded_BNL_L_LSC_Small_NaI_LY_coinc_2/figures/SDataR_Loaded_BNL_L_LSC_Small_NaI_LY_coinc_2_scale_factors.txt")

settingsFilePath = Path("/home/nick/PhD/KDK+/Large_LSC_testing/Position_test_no_outer_sleeve/2026_11_03/Position_test_no_sleeve_position_0_NaI_lead_shield/settings.xml")

savefilepath = RootDirectory / 'results'

savefilepath.mkdir(parents = True, exist_ok = True)

NaIScale = ReadInNaIScale(NaIScaleFile)

#Read in the channel labels
files = RootDirectory.rglob("*_sorted_coinc.root")

# chLabels = ReadInChannelNames(settingsFilePath)

channels = [0,1,2]

fitWindow = {'-5': [[1800,3600],[800,2100]],
             '-4': [[1600,3600],[800,2200]],
             '-3': [[1400,3200],[1000,2200]],
             '-2': [[1400,3000],[1000,2500]],
             '-1': [[1200,2800],[1400,2600]],
             '0': [[1200,2500],[1500,2500]],
             '1': [[1400,2600],[1200,2800]],
             '2': [[1000,2500],[1400,3000]],
             '3': [[1000,2200],[1500,3300]],
             '4': [[800,2200],[1500,3500]],
             '5': [[800,2000],[1500,3500]]}

Hist1DFitsCanvas = ROOT.TCanvas("Hist_1D_Position_fits", "Position Histograms",25600,1600)
ncols = 11
Hist1DFitsCanvas.Divide(ncols,2)

pData = []
for i,file in enumerate(files):
    pData.append(positionData(file,channels,ncols,fitWindow,savefilepath))  
    
    pData[-1].plot1DHist(Hist1DFitsCanvas)
    
    
Hist1DFitsCanvas.SaveAs(str(savefilepath / 'Position_1D_hist_fits.pdf'))
Hist1DFitsCanvas.SaveAs(str(savefilepath / 'Position_1D_hist_fits.png'))
