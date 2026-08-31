import ROOT
import numpy as np

def truncated_gaussian_mean_std(mu, sigma, a, b):
    """
    Mean and std of a Gaussian N(mu, sigma) truncated to [a, b].
    a, b are the raw truncation bounds (same units as mu, sigma).
    """
    fa = ROOT.TMath.Gaus(a, mu, sigma, True)   # density at lower bound
    fb = ROOT.TMath.Gaus(b, mu, sigma, True)   # density at upper bound

    cdf_a = ROOT.Math.normal_cdf(a, sigma, mu)  # note: (x, sigma, mean) order
    cdf_b = ROOT.Math.normal_cdf(b, sigma, mu)
    Z = cdf_b - cdf_a                           # retained probability mass

    mean_trunc = mu + sigma**2 * (fa - fb) / Z

    alpha = (a - mu) / sigma
    beta  = (b - mu) / sigma
    phi_alpha = ROOT.TMath.Gaus(alpha, 0, 1, True)  # standard normal density
    phi_beta  = ROOT.TMath.Gaus(beta, 0, 1, True)

    correction = (alpha * phi_alpha - beta * phi_beta) / Z - ((phi_alpha - phi_beta) / Z)**2
    var_trunc = sigma**2 * (1 + correction)
    std_trunc = var_trunc**0.5

    return mean_trunc, std_trunc

def propagate_truncated_mean_error(mu, sigma, a, b, cov_matrix, h=1e-5):
    """
    cov_matrix: 2x2 covariance matrix for (mu, sigma), e.g. from TFitResultPtr
    h: step size for numerical derivative
    """
    def f(m, s):
        return truncated_gaussian_mean_std(m, s, a, b)[0]

    def sig(m,s):
        return truncated_gaussian_mean_std(m,s,a,b)[1]

    dmean_dmu    = (f(mu+h, sigma) - f(mu-h, sigma)) / (2*h)
    dmean_dsigma = (f(mu, sigma+h) - f(mu, sigma-h)) / (2*h)

    var_mu, var_sigma = cov_matrix(0,0), cov_matrix(1,1)
    cov_mu_sigma = cov_matrix(0,1)

    variance = (dmean_dmu**2 * var_mu
                + dmean_dsigma**2 * var_sigma
                + 2 * dmean_dmu * dmean_dsigma * cov_mu_sigma)
    
    
    dmean_dmu    = (sig(mu+h, sigma) - sig(mu-h, sigma)) / (2*h)
    dmean_dsigma = (sig(mu, sigma+h) - sig(mu, sigma-h)) / (2*h)

    var_mu, var_sigma = cov_matrix(0,0), cov_matrix(1,1)
    cov_mu_sigma = cov_matrix(0,1)

    Sigvariance = (dmean_dmu**2 * var_mu
                + dmean_dsigma**2 * var_sigma
                + 2 * dmean_dmu * dmean_dsigma * cov_mu_sigma)
    return variance**0.5, Sigvariance**0.5

### Define the initial parameters of the function, then adjust them through each iteration.
truncRange = [10,150]
threshold = 60
mu = 80
sigma = 15
sqrtRunNumber = 5

truncFormula = (
    f"[0]*TMath::Gaus(x,[1],[2],1)/"
    f"(ROOT::Math::normal_cdf({threshold},[2],[1]) - ROOT::Math::normal_cdf({truncRange[0]},[2],[1]))"
)
func = ROOT.TF1(f"truncGaus", truncFormula, threshold, truncRange[1])

GausFit = ROOT.TF1("GaussianFit","gaus",truncRange[0],truncRange[1])
GausFit.SetLineColor(ROOT.kRed)

func.SetParameters(1.0,mu,sigma) #Amplitude, Mean, Sigma
func.SetLineColor(ROOT.kGreen+2)

canvas = ROOT.TCanvas("c1","Decreasing Mean Cumulative Data", 3200,800)
canvas.Divide(sqrtRunNumber,sqrtRunNumber)

noDecCumCanvas = ROOT.TCanvas("c2", "No Decrease Cummulative", 3200, 800)
noDecCumCanvas.Divide(sqrtRunNumber,sqrtRunNumber)

GaussianData = ROOT.TH1F("Sample_Gaussian","Sampled Gaussian",100,truncRange[0],truncRange[1]) #Define the Gaussian Data and then keep filling it iteratively. 
GaussianData.Reset() 

noDecGaussianData = ROOT.TH1F("Sample_Gaussian_no_decrease","Sampled Gaussian No Decrease",100,truncRange[0],truncRange[1]) #Define the Gaussian Data and then keep filling it iteratively. 
noDecGaussianData.Reset() 

snapshots = []  
genFuncs = [] 
RMeans,GMeans,TMeans = [],[],[]
RStd,GStd,TStd = [],[],[]

RMeansErr,GMeansErr,TMeansErr = [],[],[]
RStdErr,GStdErr,TStdErr = [],[],[]

iterMeans = []
TmeanLines,GmeanLines = [],[]


noDecsnapshots = []  
noDecgenFuncs = [] 
noDecRMeans,noDecGMeans,noDecTMeans = [],[],[]
noDecRStd,noDecGStd,noDecTStd = [],[],[]

noDecRMeansErr,noDecGMeansErr,noDecTMeansErr = [],[],[]
noDecRStdErr,noDecGStdErr,noDecTStdErr = [],[],[]

noDeciterMeans = []
noDecTmeanLines,noDecGmeanLines = [],[]

for i in range(sqrtRunNumber**2):
    pad = canvas.cd(i+1)
    try:
        GaussData = ROOT.TF1(f"Gaussian_Func_{i}","gaus",truncRange[0],truncRange[1])
        mean = mu - i
        iterMeans.append(mean)
        GaussData .SetParameters(1.0,mean,sigma) #Amplitude, Mean, Sigma
        genFuncs.append(GaussData)
        

        nSamples = 10000

        for j in range(nSamples):
            var = GaussData.GetRandom()
            if var > threshold:
                GaussianData.Fill(var)
            
        snapshot = GaussianData.Clone(f"snapshot_{i}")
        snapshot.SetDirectory(0)  # detach from any directory, avoid double-registration
        snapshots.append(snapshot) 
        snapshot.Draw("same")

        # GaussianData.Draw("same")
        GausResults = snapshot.Fit(GausFit,"RSQ")
        TruncResults = snapshot.Fit(func,"RSQ+")
        
        RMeans.append(snapshot.GetMean())
        RStd.append(snapshot.GetStdDev())
        
        GMeans.append(GausFit.GetParameter(1))
        GMeansErr.append(GausFit.GetParError(1))
        GStd.append(GausFit.GetParameter(2))
        GStdErr.append(GausFit.GetParError(2))
        
        TMeans.append(func.GetParameter(1))
        TMeansErr.append(func.GetParError(1))
        TStd.append(func.GetParameter(2))
        TStdErr.append(func.GetParError(2))
        
        
        ymax = snapshot.GetMaximum() * 1.05
        
        TmeanLines.append(ROOT.TLine(TMeans[-1],0,TMeans[-1],ymax))
        TmeanLines[-1].SetLineColor(ROOT.kGreen)
        TmeanLines[-1].SetLineWidth(2)
        TmeanLines[-1].SetLineStyle(2)
        TmeanLines[-1].Draw("SAME")
        
        GmeanLines.append(ROOT.TLine(GMeans[-1],0,GMeans[-1],ymax))
        GmeanLines[-1].SetLineColor(ROOT.kRed)
        GmeanLines[-1].SetLineWidth(2)
        GmeanLines[-1].SetLineStyle(2)
        GmeanLines[-1].Draw("SAME")
        

        canvas.Update()
    except Exception as e:
        print(f"Failed on iteration {i}: {e}")
    
    #############################################################################################
    ########### Make a cummulative distribution that doesn't have a decreasing mean.    #########
     #############################################################################################
    pad = noDecCumCanvas.cd(i+1)
    try:
        GaussData = ROOT.TF1(f"Gaussian_Func_{i}_no_decrease","gaus",truncRange[0],truncRange[1])
        mean = mu 
        noDeciterMeans.append(mean)
        GaussData.SetParameters(1.0,mean,sigma) #Amplitude, Mean, Sigma
        noDecgenFuncs.append(GaussData)
        

        nSamples = 10000

        for j in range(nSamples):
            var = GaussData.GetRandom()
            if var > threshold:
                noDecGaussianData.Fill(var)
            
        snapshot = noDecGaussianData.Clone(f"snapshot_{i}")
        snapshot.SetDirectory(0)  # detach from any directory, avoid double-registration
        noDecsnapshots.append(snapshot) 
        snapshot.Draw("same")

        # GaussianData.Draw("same")
        GausResults = snapshot.Fit(GausFit,"RSQ")
        TruncResults = snapshot.Fit(func,"RSQ+")
        cov = TruncResults.GetCovarianceMatrix() 
        
        noDecRMeans.append(snapshot.GetMean())
        noDecRStd.append(snapshot.GetStdDev())
        
        noDecGMeans.append(GausFit.GetParameter(1))
        noDecGMeansErr.append(GausFit.GetParError(1))
        noDecGStd.append(GausFit.GetParameter(2))
        noDecGStdErr.append(GausFit.GetParError(2))
        
        Tmu = func.GetParameter(1)
        TSigma = func.GetParameter(2)
        TMean, TStdTrunc = truncated_gaussian_mean_std(Tmu, TSigma, threshold, truncRange[1])
        TMeanErr,TStdTruncErr = propagate_truncated_mean_error(Tmu, TSigma, threshold, truncRange[1], cov)

        '''
        truncFormula = (
        f"[0]*TMath::Gaus(x,[1],[2],1)/"
        f"(ROOT::Math::normal_cdf({truncRange[1]},[2],[1]) - ROOT::Math::normal_cdf({truncRange[0]},[2],[1]))"
)
        '''
        
        
        noDecTMeans.append(func.GetParameter(1))
        # noDecTMeans.append(TMean)
        noDecTMeansErr.append(func.GetParError(1))
        # noDecTMeansErr.append(TMeanErr)
        noDecTStd.append(func.GetParameter(2))
        # noDecTStd.append(TStdTrunc)
        noDecTStdErr.append(func.GetParError(2))
        # noDecTStdErr.append(TStdTruncErr)
        
        
        ymax = snapshot.GetMaximum() * 1.05
        
        noDecTmeanLines.append(ROOT.TLine(noDecTMeans[-1],0,noDecTMeans[-1],ymax))
        noDecTmeanLines[-1].SetLineColor(ROOT.kGreen)
        noDecTmeanLines[-1].SetLineWidth(2)
        noDecTmeanLines[-1].SetLineStyle(2)
        noDecTmeanLines[-1].Draw("SAME")
        
        noDecGmeanLines.append(ROOT.TLine(noDecGMeans[-1],0,noDecGMeans[-1],ymax))
        noDecGmeanLines[-1].SetLineColor(ROOT.kRed)
        noDecGmeanLines[-1].SetLineWidth(2)
        noDecGmeanLines[-1].SetLineStyle(2)
        noDecGmeanLines[-1].Draw("SAME")
        

        noDecCumCanvas.Update()
    except Exception as e:
        print(f"Failed on iteration {i}: {e}")

# canvas.Draw()

print("Decreasing Mean Cummulative Data")
for i in range(sqrtRunNumber**2):
    print(f"Iteration Number: {i}\t Root Mean: {RMeans[i]}\t Gaus Fit Mean: {GMeans[i]}\t Truncated Fit Mean: {TMeans[i]}\t Iteration Mean: {iterMeans[i]}")
    
print()
for i in range(sqrtRunNumber**2):
    print(f"Iteration Number: {i}\t Root STD: {RStd[i]}\t Gaus Fit STD: {GStd[i]}\t Truncated Fit: {TStd[i]}")
    
print("Non-Decreasing Mean Cummulative Data")
for i in range(sqrtRunNumber**2):
    print(f"Iteration Number: {i}\t Root Mean: {noDecRMeans[i]}\t Gaus Fit Mean: {noDecGMeans[i]}\t Truncated Fit Mean: {noDecTMeans[i]}\t Iteration Mean: {noDeciterMeans[i]}")
    
print()
for i in range(sqrtRunNumber**2):
    print(f"Iteration Number: {i}\t Root STD: {noDecRStd[i]}\t Gaus Fit STD: {noDecGStd[i]}\t Truncated Fit: {noDecTStd[i]}")
    

    
    
    
resCanvas = ROOT.TCanvas("STD_Mean_Res_results", "STD Mean Resolution Results", 1600,800)

resCanvas.Divide(3,2)

index = np.array([i for i in range(len(TMeans))], dtype  = 'float64')
indexErr = np.zeros(len(TMeans),dtype='float64')

TMeanArr = np.array(TMeans,dtype='float64')
TMeanErrArr = np.array(TMeansErr,dtype='float64')

TStdArr = np.array(TStd,dtype='float64')
TStdErrArr = np.array(TStdErr,dtype='float64')

TRes = np.array([s/m for s,m in zip(TStdArr,TMeanArr)],dtype='float64')
TResErr = np.array([np.sqrt((s/m)**2 * ((sErr/s)**2 + (mErr/m)**2)) for s,sErr,m,mErr in zip(TStdArr,TStdErrArr,TMeanArr,TMeanErrArr)], dtype='float64')


GMeanArr = np.array(GMeans,dtype='float64')
GMeanErrArr = np.array(GMeansErr,dtype='float64')

GStdArr = np.array(GStd,dtype='float64')
GStdErrArr = np.array(GStdErr,dtype='float64')

GRes = np.array([s/m for s,m in zip(GStdArr,GMeanArr)],dtype='float64')
GResErr = np.array([np.sqrt((s/m)**2 * ((sErr/s)**2 + (mErr/m)**2)) for s,sErr,m,mErr in zip(GStdArr,GStdErrArr,GMeanArr,GMeanErrArr)], dtype='float64')


pad = resCanvas.cd(1)
GMeanDist = ROOT.TGraphErrors(len(GMeans),index,GMeanArr,indexErr,GMeanErrArr)
pad.SetLeftMargin(0.25)
GMeanDist.SetMarkerStyle(20)
GMeanDist.SetMarkerSize(1.0)
GMeanDist.SetMarkerColor(ROOT.kBlue)
GMeanDist.GetXaxis().SetTitle("Day Number")
GMeanDist.GetYaxis().SetTitle("Gaussian Mean")
GMeanDist.SetTitle("Gaussian Mean")

GMeanDist.Draw("AP")
resCanvas.Update()

pad = resCanvas.cd(2)
GStdDist = ROOT.TGraphErrors(len(GStd),index,GStdArr,indexErr,GStdErrArr)
pad.SetLeftMargin(0.25)
GStdDist.SetMarkerStyle(20)
GStdDist.SetMarkerSize(1.0)
GStdDist.SetMarkerColor(ROOT.kGreen)
GStdDist.GetXaxis().SetTitle("Day Number")
GStdDist.GetYaxis().SetTitle("Gaussian STD")
GStdDist.SetTitle("Gaussian Std")


GStdDist.Draw("AP")
resCanvas.Update()

pad = resCanvas.cd(3)
GResDist = ROOT.TGraphErrors(len(GMeans),index,GRes,indexErr,GResErr)
pad.SetLeftMargin(0.25)
GResDist.SetMarkerStyle(20)
GResDist.SetMarkerSize(1.0)
GResDist.SetMarkerColor(ROOT.kRed)
GResDist.GetXaxis().SetTitle("Day Number")
GResDist.GetYaxis().SetTitle("Gaussian resolution")
GResDist.SetTitle("Gaussian Resolution")

GResDist.Draw("AP")
resCanvas.Update()




pad = resCanvas.cd(4)
TMeanDist = ROOT.TGraphErrors(len(TMeans),index,TMeanArr,indexErr,TMeanErrArr)
pad.SetLeftMargin(0.25)
TMeanDist.SetMarkerStyle(20)
TMeanDist.SetMarkerSize(1.0)
TMeanDist.SetMarkerColor(ROOT.kBlue)
TMeanDist.GetXaxis().SetTitle("Day Number")
TMeanDist.GetYaxis().SetTitle("Truncated Gaussian Mean")
TMeanDist.SetTitle("Truncated Gaussian Mean")

TMeanDist.Draw("AP")
resCanvas.Update()

pad = resCanvas.cd(5)
TStdDist = ROOT.TGraphErrors(len(GStd),index,TStdArr,indexErr,TStdErrArr)
pad.SetLeftMargin(0.25)
TStdDist.SetMarkerStyle(20)
TStdDist.SetMarkerSize(1.0)
TStdDist.SetMarkerColor(ROOT.kGreen)
TStdDist.GetXaxis().SetTitle("Day Number")
TStdDist.GetYaxis().SetTitle("Truncated Gaussian STD")
TStdDist.SetTitle("Truncated Gaussian Std")

TStdDist.Draw("AP")
resCanvas.Update()

pad = resCanvas.cd(6)
TResDist = ROOT.TGraphErrors(len(GMeans),index,TRes,indexErr,TResErr)
pad.SetLeftMargin(0.25)
TResDist.SetMarkerStyle(20)
TResDist.SetMarkerSize(1.0)
TResDist.SetMarkerColor(ROOT.kRed)
TResDist.GetXaxis().SetTitle("Day Number")
TResDist.GetYaxis().SetTitle("Truncated Gaussian Resolution")
TResDist.SetTitle("Truncated Gaussian Resolution")

TResDist.Draw("AP")
resCanvas.Update()



###########################################################################################################
#################   Plot the non-decreasing cumulative Results  ###########################################
###########################################################################################################



noDecresCanvas = ROOT.TCanvas("STD_Mean_Res_results_no_dec", "STD Mean Resolution No Decreasing Mean Results", 1600,800)

noDecresCanvas.Divide(3,2)

index = np.array([i for i in range(len(TMeans))], dtype  = 'float64')
indexErr = np.zeros(len(TMeans),dtype='float64')

noDecTMeanArr = np.array(noDecTMeans,dtype='float64')
noDecTMeanErrArr = np.array(noDecTMeansErr,dtype='float64')

noDecTStdArr = np.array(noDecTStd,dtype='float64')
noDecTStdErrArr = np.array(noDecTStdErr,dtype='float64')

noDecTRes = np.array([s/m for s,m in zip(noDecTStdArr,noDecTMeanArr)],dtype='float64')
noDecTResErr = np.array([np.sqrt((s/m)**2 * ((sErr/s)**2 + (mErr/m)**2)) for s,sErr,m,mErr in zip(noDecTStdArr,noDecTStdErrArr,noDecTMeanArr,noDecTMeanErrArr)], dtype='float64')


noDecGMeanArr = np.array(noDecGMeans,dtype='float64')
noDecGMeanErrArr = np.array(noDecGMeansErr,dtype='float64')

noDecGStdArr = np.array(noDecGStd,dtype='float64')
noDecGStdErrArr = np.array(noDecGStdErr,dtype='float64')

noDecGRes = np.array([s/m for s,m in zip(noDecGStdArr,noDecGMeanArr)],dtype='float64')
noDecGResErr = np.array([np.sqrt((s/m)**2 * ((sErr/s)**2 + (mErr/m)**2)) for s,sErr,m,mErr in zip(noDecGStdArr,noDecGStdErrArr,noDecGMeanArr,noDecGMeanErrArr)], dtype='float64')


pad = noDecresCanvas.cd(1)
noDecGMeanDist = ROOT.TGraphErrors(len(GMeans),index,noDecGMeanArr,indexErr,noDecGMeanErrArr)
pad.SetLeftMargin(0.25)
noDecGMeanDist.SetMarkerStyle(20)
noDecGMeanDist.SetMarkerSize(1.0)
noDecGMeanDist.SetMarkerColor(ROOT.kBlue)
noDecGMeanDist.GetXaxis().SetTitle("Day Number")
noDecGMeanDist.GetYaxis().SetTitle("Gaussian Mean")
noDecGMeanDist.SetTitle("Gaussian Mean")

noDecGMeanDist.Draw("AP")
noDecresCanvas.Update()

pad = noDecresCanvas.cd(2)
noDecGStdDist = ROOT.TGraphErrors(len(GStd),index,noDecGStdArr,indexErr,noDecGStdErrArr)
pad.SetLeftMargin(0.25)
noDecGStdDist.SetMarkerStyle(20)
noDecGStdDist.SetMarkerSize(1.0)
noDecGStdDist.SetMarkerColor(ROOT.kGreen)
noDecGStdDist.GetXaxis().SetTitle("Day Number")
noDecGStdDist.GetYaxis().SetTitle("Gaussian STD")
noDecGStdDist.SetTitle("Gaussian Std")


noDecGStdDist.Draw("AP")
noDecresCanvas.Update()

pad = noDecresCanvas.cd(3)
noDecGResDist = ROOT.TGraphErrors(len(GMeans),index,noDecGRes,indexErr,noDecGResErr)
pad.SetLeftMargin(0.25)
noDecGResDist.SetMarkerStyle(20)
noDecGResDist.SetMarkerSize(1.0)
noDecGResDist.SetMarkerColor(ROOT.kRed)
noDecGResDist.GetXaxis().SetTitle("Day Number")
noDecGResDist.GetYaxis().SetTitle("Gaussian resolution")
noDecGResDist.SetTitle("Gaussian Resolution")

noDecGResDist.Draw("AP")
noDecresCanvas.Update()




pad = noDecresCanvas.cd(4)
noDecTMeanDist = ROOT.TGraphErrors(len(TMeans),index,noDecTMeanArr,indexErr,noDecTMeanErrArr)
pad.SetLeftMargin(0.25)
noDecTMeanDist.SetMarkerStyle(20)
noDecTMeanDist.SetMarkerSize(1.0)
noDecTMeanDist.SetMarkerColor(ROOT.kBlue)
noDecTMeanDist.GetXaxis().SetTitle("Day Number")
noDecTMeanDist.GetYaxis().SetTitle("Truncated Gaussian Mean")
noDecTMeanDist.SetTitle("Truncated Gaussian Mean")

noDecTMeanDist.Draw("AP")
noDecresCanvas.Update()

pad = noDecresCanvas.cd(5)
noDecTStdDist = ROOT.TGraphErrors(len(GStd),index,noDecTStdArr,indexErr,noDecTStdErrArr)
pad.SetLeftMargin(0.25)
noDecTStdDist.SetMarkerStyle(20)
noDecTStdDist.SetMarkerSize(1.0)
noDecTStdDist.SetMarkerColor(ROOT.kGreen)
noDecTStdDist.GetXaxis().SetTitle("Day Number")
noDecTStdDist.GetYaxis().SetTitle("Truncated Gaussian STD")
noDecTStdDist.SetTitle("Truncated Gaussian Std")

noDecTStdDist.Draw("AP")
noDecresCanvas.Update()

pad = noDecresCanvas.cd(6)
noDecTResDist = ROOT.TGraphErrors(len(GMeans),index,noDecTRes,indexErr,noDecTResErr)
pad.SetLeftMargin(0.25)
noDecTResDist.SetMarkerStyle(20)
noDecTResDist.SetMarkerSize(1.0)
noDecTResDist.SetMarkerColor(ROOT.kRed)
noDecTResDist.GetXaxis().SetTitle("Day Number")
noDecTResDist.GetYaxis().SetTitle("Truncated Gaussian Resolution")
noDecTResDist.SetTitle("Truncated Gaussian Resolution")

noDecTResDist.Draw("AP")
noDecresCanvas.Update()