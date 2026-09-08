import numpy as np
import ROOT
from pathlib import Path
import pandas as pd

class sourceData:
    def __init__(self,source,file,fitWindows,nbins):
        self.source = source
        self.file = file
        self.fitWindows = fitWindows
        self.fitFunc = None
       

        self.NaI1Data = ROOT.TH1D(f"NaI_1_Hist_data_{source}",f"NaI 1 Hist Data {source}",nbins,0,4000)
        self.NaI2Data = ROOT.TH1D(f"NaI_2_Hist_data_{source}",f"NaI 2 Hist Data {source}",nbins,0,4000)
        self.NaI3Data = ROOT.TH1D(f"NaI_3_Hist_data_{source}",f"NaI 3 Hist Data {source}",nbins,0,4000)
        self.NaI4Data = ROOT.TH1D(f"NaI_4_Hist_data_{source}",f"NaI 4 Hist Data {source}",nbins,0,4000)
        
        self.ChDict = {'8': [self.NaI1Data,self.fitWindows[0]],
                       '10': [self.NaI2Data,self.fitWindows[1]],
                       '12': [self.NaI3Data,self.fitWindows[2]],
                       '14': [self.NaI4Data,self.fitWindows[3]]}
        
    def ReadInData(self):
        #Reads in the data from the file and then saves it to the appropriate histogram
        df = pd.read_csv(self.file,sep = ";",skiprows = 1, header = None,
                         usecols = [1,3],names = ["ch", "Energy"], dtype = {"ch": "int32","Energy": "float64"},engine = "pyarrow")
        
        for ch,hist in self.ChDict.items():
            chValue = int(ch)
            Energy = df.loc[df["ch"] == chValue, "Energy"].to_numpy()
            hist[0].FillN(len(Energy),Energy,np.ones(len(Energy)))

    def FitData(self):    
        
         
        
        self.DataCanvas = ROOT.TCanvas(f"{self.source}_canvas",f"{self.source} Canvas", 3200,1600)
        self.DataCanvas.Divide(2,2)
        
            
        padInd = 1  
        for key,value in self.ChDict.items():
            
            hist = value[0]
            fitWindow = value[1]
            
            if self.source == "Co-60":
                # fitFormula = ("[0]*TMath::Gaus(x,[1],[2],1) + [3]*TMath::Gaus(x,[4],[5],1) + [6]*TMath::Expo(x)")
                self.fitFunc =  ROOT.TF1(f"{self.source}_{key}_double_gaus","gaus(0) + gaus(3) + expo(6)", fitWindow[0], fitWindow[1])
                ampGuess = 5000
                sigGuess = fitWindow[1]-fitWindow[0]
                muGuess = (fitWindow[1]+fitWindow[0])/3
                muGuess2 = 2*(fitWindow[1]+fitWindow[0])/3
                
                expAmpGuess = 5
                tau = 1e-10
                self.fitFunc.SetParameters(ampGuess,muGuess,sigGuess,ampGuess,muGuess2,sigGuess,expAmpGuess,tau)
                self.fitFunc.SetLineColor(ROOT.kRed)
                self.fitFunc.SetLineWidth(5)
                
            elif self.source == "Na-22":
                self.fitFunc1 =  ROOT.TF1(f"{self.source}_{key}_double_gaus_1","gaus(0) + expo(3)", fitWindow[0], fitWindow[1])
                self.fitFunc2 =  ROOT.TF1(f"{self.source}_{key}_double_gaus_2","gaus(0) + expo(3)", fitWindow[2], fitWindow[3])
                self.fitFunc3 =  ROOT.TF1(f"{self.source}_{key}_double_gaus_3","gaus(0) + expo(3)", fitWindow[4], fitWindow[5])
                ampGuess = 1000
                expAmpGuess = 5
                sigGuess1 = fitWindow[1]-fitWindow[0]
                muGuess1 = (fitWindow[1]+fitWindow[0])/2
                
                sigGuess2 = fitWindow[3]-fitWindow[2]
                muGuess2 = (fitWindow[3]+fitWindow[2])/2
                
                sigGuess3 = fitWindow[5]-fitWindow[4]
                muGuess3 = (fitWindow[5]+fitWindow[4])/2
                tau = 1e-10
                self.fitFunc1.SetParameters(ampGuess,muGuess1,sigGuess1,expAmpGuess,tau)
                self.fitFunc2.SetParameters(ampGuess,muGuess2,sigGuess2,expAmpGuess,tau)
                self.fitFunc3.SetParameters(ampGuess,muGuess3,sigGuess3,expAmpGuess,tau)
                
                self.fitFunc1.SetLineColor(ROOT.kRed)
                self.fitFunc2.SetLineColor(ROOT.kRed)
                self.fitFunc3.SetLineColor(ROOT.kRed)
                
                self.fitFunc1.SetLineWidth(5)
                self.fitFunc2.SetLineWidth(5)
                self.fitFunc3.SetLineWidth(5)
                
            else:
                self.fitFunc =  ROOT.TF1(f"{self.source}_{key}_double_gaus","gaus(0) + expo(3)", fitWindow[0], fitWindow[1])
                ampGuess = 5000
                sigGuess = fitWindow[1]-fitWindow[0]
                muGuess = (fitWindow[1]+fitWindow[0])/2
                expAmpGuess = 5
                tau = 1e-10
                self.fitFunc.SetParameters(ampGuess,muGuess,sigGuess,expAmpGuess,tau)
                self.fitFunc.SetLineColor(ROOT.kRed)
                self.fitFunc.SetLineWidth(5)
            
            
            
            
            pad = self.DataCanvas.cd(padInd)
            pad.SetLogy(1)
            
            hist.Draw("Same")
            if self.fitFunc is not None:
                hist.Fit(self.fitFunc,"R")
                self.fitFunc.Draw("SAME")
            else:
                hist.Fit(self.fitFunc1,"R")
                hist.Fit(self.fitFunc2,"R")
                hist.Fit(self.fitFunc3,"R")
                
                self.fitFunc1.Draw("SAME")
                self.fitFunc2.Draw("SAME")
                self.fitFunc3.Draw("SAME")
            
            self.DataCanvas.Update()

                
            
            padInd +=1
        
                
                
def ReadInInitFile(file,nbins):
    
    
    with open(file) as f:
        
        lines = f.readlines()
        
        saveFilePath = lines[0].split("\t")[1]
        
        data = []
        
        for i in range(5):
            source = lines[i+1].split("\t")[0]
            sourceFile = lines[i+1].split("\n")[0].split("\t")[1]
            fitWindowsStr = lines[i+6].split("\n")[0].split("\t")[1:]
            
            fitWindowsInt = [[int(x) for x in item.split(",")] for item in fitWindowsStr]
            
            data.append(sourceData(source,sourceFile,fitWindowsInt,nbins))
            
    return data
        
    
                
            
initFile = Path("/home/nick/PhD/KDK+/Annulus_Compton_scatter_V1/NaI_Annulus_calibration/CFD_lower_HV_Inside_Source.txt")   

nbins = 100

sourceData = ReadInInitFile(initFile,nbins)

for data in sourceData:
    print(f"Reading in {data.source} data:")
    data.ReadInData()
    data.FitData()
 
        

        
        