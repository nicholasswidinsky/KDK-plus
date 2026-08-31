import ROOT
import numpy as np


class Hist2DData:
    def __init__(self,threshold,plotRange,rotatedPlotRange,fitRange,totE,minE,maxE,xRes,yRes):
        
        self.totE = totE
        self.minE = minE
        self.maxE = maxE
        
        self.xRes = xRes
        self.yRes = yRes
        
        self.xValues = []
        self.yValues = []
        
        self.s0 = 5.0 #20 keV, sigma at 0 energy
        self.s1 = 0.12 #keV per keV, linear scale
        self.s2 = 8.0e-6 #keV per keV^2, quadratic scale factor
        self.sSqrt = 3 #keV per sqrt(keV), sqrt scale factor
        
        self.nbins = 100
        
        self.threshold = threshold
        self.plotRange = plotRange
        self.rotatedPlotRange = rotatedPlotRange
        self.fitRange = fitRange
        
        self.hist2DData = []
        self.RotatedHist2DData = []
        self.SummedInt = []

        truncFormula = (f"[0]*TMath::Gaus(x,[1],[2],1)/"
                            f"(ROOT::Math::normal_cdf({self.threshold},[2],[1]) - ROOT::Math::normal_cdf({self.plotRange[1]},[2],[1]))")
        
        self.truncFunc = ROOT.TF1(f"TruncGaus", truncFormula,self.threshold,self.plotRange[1])
        self.truncFunc.SetLineColor(ROOT.kGreen)
        
        self.GausFunc = ROOT.TF1(f"GaussianFunction", "gaus", self.plotRange[0],self.plotRange[1])
        self.GausFunc.SetLineColor(ROOT.kRed)
        
    def GenData(self,nSamples,decreaseStep,iterations,DecreaseX = False):
        # This function will be called every iteration to generate the data for that iteration. 
        
        # coincHist = ROOT.TH2D(f"2D_histogram", f"2D Histogram", self.nbins, self.plotRange[0],self.plotRange[1],self.nbins, self.plotRange[0],self.plotRange[1])
                    
        for i in range(iterations):
            
            self.xValues.append([])
            self.yValues.append([])
            #I think that if I use this version of the coincHist object I should be able to just plot the 2D hist for every day, rather than the cummulative data. 
            
            coincHist = ROOT.TH2D(f"2D_histogram_iteration_{i}", f"2D Histogram iteration {i}", self.nbins, self.plotRange[0],self.plotRange[1],self.nbins, self.plotRange[0],self.plotRange[1])
            
            # if DecreaseX:
            #     #If I want to decrease the energy of the LSC I just decrease what the total energy allowed it to grab from the uniform dist.
            #     #Since I am always sampling the x axis first (LSC) I don't need to worry about changing the value of the max E.
            #     self.maxE = self.maxE - decreaseStep
            
            for j in range(nSamples):

                # xValue= ROOT.gRandom.Gaus(500,100)
                xValue= ROOT.gRandom.Uniform(self.minE,self.maxE)
                yValue = self.totE - xValue
                
                # xRes = self.s0 + np.sqrt(xValue) * self.sSqrt
                xRes = self.s0 + self.s1*xValue
                # xRes = self.xRes
                
                #Apply a Gaussian Smear based on the resolution of the detectors. 
                if DecreaseX:
                    #Set the LSC to lose a set percent energy for every iteration. This is still a rough estimate and doesn't fully respresent the system. 
                    redXvalue = (xValue * (1-(decreaseStep * i)))
                    # redXvalue = xValue #No decrease in the xValue here for testing purposes. 
                    if redXvalue < self.minE:
                        redXvalue = self.minE
                    xVar = redXvalue  + ROOT.gRandom.Gaus(0,xRes)

                else:
                    xVar = xValue + ROOT.gRandom.Gaus(0,xRes)
                    
                yVar = yValue + ROOT.gRandom.Gaus(0,self.yRes)
                if xVar + yVar > self.threshold:
                    self.xValues[i].append(xVar)
                    self.yValues[i].append(yVar)
                    coincHist.Fill(xVar,yVar)

            iterData = coincHist.Clone(f"iteration_{i}")
            iterData.SetDirectory(0)
            self.hist2DData.append(iterData) #Save the cumulative data that I just generated. 
             
    def rotateHist(self):
        if len(self.hist2DData) == 0:
            print('Data has not been generated')
        else:
            #Use this definition if I am plotting cumulative. 
            rotatedcoincHist = ROOT.TH2D(f"2D_histogram_rotated", f"2D Histogram Rotated", self.nbins, self.rotatedPlotRange[0],self.rotatedPlotRange[1],self.nbins, self.plotRange[0],self.plotRange[1])
            
            SummedInt = ROOT.TH1D(f"Summed_Int_cummulative", "Summed Integral cummulative", self.nbins, self.plotRange[0],self.plotRange[1])
            
            i = 0
            for xData,yData in zip(self.xValues,self.yValues):
                # rotatedcoincHist = ROOT.TH2D(f"2D_histogram_rotated_iteration_{i}", f"2D Histogram Rotated iteration {i}", self.nbins, self.rotatedPlotRange[0],self.rotatedPlotRange[1],self.nbins, self.plotRange[0],self.plotRange[1])
                
                # SummedInt = ROOT.TH1D(f"Summed_Int_noncummulative_iteration_{i}", f"Summed Integral non-cummulative iteration {i}", self.nbins, self.plotRange[0],self.plotRange[1])
                
                i +=1 
                for x,y in zip(xData,yData):
                    sumData = x + y
                    difData = x-y
                    
                    SummedInt.Fill(sumData)
                    rotatedcoincHist.Fill(difData,sumData)
                    
                    
                iterData = rotatedcoincHist.Clone(f"rot_iteration_{i}")
                iterData.SetDirectory(0)
                self.RotatedHist2DData.append(iterData)
                
                iterSumData = SummedInt.Clone(f"Summed_int_iteration_{i}")
                iterSumData.SetDirectory(0)
                self.SummedInt.append(iterSumData)
                
                
    def Plot2DHist(self,nCols):
        
        self.hist2DCanvas = ROOT.TCanvas("Hist_2D_Canvas","Hist 2D Canvas", 3200, 1600)
        self.hist2DCanvas.Divide(nCols,nCols)
        
        for i in range(nCols**2):
            
            pad = self.hist2DCanvas.cd(i+1)
            
            self.hist2DData[i].Draw("SAME")
            
            self.hist2DCanvas.Update()
            
    def PlotRotated2DHist(self,nCols):
        self.Rothist2DCanvas = ROOT.TCanvas("Rotated_Hist_2D_Canvas","Rotated Hist 2D Canvas", 3200, 1600)
        self.Rothist2DCanvas.Divide(nCols,nCols)
        
        for i in range(nCols**2):
            
            pad = self.Rothist2DCanvas.cd(i+1)
            
            self.RotatedHist2DData[i].Draw("SAME")
            
            self.Rothist2DCanvas.Update()
            
    def PlotSummedInt(self,nCols):
        self.SummedIntCanvas = ROOT.TCanvas("Summed_int_hist", "Summed Integral Histogram", 3200,1600)
        self.SummedIntCanvas.Divide(nCols,nCols)
        
        dGausfunc = ROOT.TF1(f"Double_gaus","gaus(0) + gaus(3)",self.fitRange[0],self.fitRange[1])
        dGausfunc.SetLineColor(ROOT.kRed)
        dGausfunc.SetParLimits(0,0.0,1e4)
        dGausfunc.SetParLimits(1,self.fitRange[0],self.fitRange[1])
        dGausfunc.SetParLimits(2,0.0,1e4)
        dGausfunc.SetParLimits(3,0.0,1e4)
        dGausfunc.SetParLimits(4,self.fitRange[0],self.fitRange[1])
        dGausfunc.SetParLimits(5,0.0,1e4)
        
        Gausfunc = ROOT.TF1("gaus_func", "gaus(0)", self.fitRange[0],self.fitRange[1])
        Gausfunc.SetLineColor(ROOT.kGreen)
        Gausfunc.SetParLimits(0,0.0,1e10)
        Gausfunc.SetParLimits(1,self.fitRange[0],self.fitRange[1])
        Gausfunc.SetParLimits(2,0.0,1e10)
        
        
        lorentzFormula = "[0]*TMath::CauchyDist(x,[1],[2])"
        lorentzFit = ROOT.TF1("lorentz_fit",lorentzFormula,self.plotRange[0],self.plotRange[1])
        lorentzFit.SetLineColor(ROOT.kMagenta)

        
        for i in range(nCols**2):
            pad = self.SummedIntCanvas.cd(i+1)
            
            amp = self.SummedInt[i].Integral('Width')
            mu = self.SummedInt[i].GetMean()
            sig = self.SummedInt[i].GetStdDev()
            
            print(f"Amplitude Guess: {amp}")
            print(f"Mu Guess: {mu}")
            print(f"Sigma Guess: {sig}")
            
            dGausfunc.SetParameters(amp*0.05,mu,sig,amp*0.05,mu,sig)
            Gausfunc.SetParameters(amp,mu,sig)
            lorentzFit.SetParameters(amp,mu,sig*2)
            
            self.SummedInt[i].Draw("SAME")
            self.SummedInt[i].Fit(dGausfunc,"R")
            self.SummedInt[i].Fit(Gausfunc,"R+")
            self.SummedInt[i].Fit(lorentzFit,"R+")
            self.SummedIntCanvas.Update()
       
threshold = 0
minE = 75
maxE = 550
totE = 662   
plotRange = [0,900]
rotatedPlotRange = [-600,600]
fitRange = [550,750]
xRes = 50 #LSC resolution ~ 15%
yRes = 5 #NaI resolution ~ 7%
nSamples = 10000
decreaseStep = 0.000 #LSC Looses this percent of LY each iteration.
sqrtIterations = 2
iterations = sqrtIterations **2

data = Hist2DData(threshold,plotRange,rotatedPlotRange,fitRange,totE,minE,maxE,xRes,yRes)   
data.GenData(nSamples,decreaseStep,iterations,DecreaseX = True)  
data.Plot2DHist(sqrtIterations)
data.rotateHist()
data.PlotRotated2DHist(sqrtIterations)
data.PlotSummedInt(sqrtIterations)
