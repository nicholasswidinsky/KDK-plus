import numpy as np
from pathlib import Path


openedfiles = [] #List of all the files that have been opened. This allows me to write the header the first time, then just append to the file every other time. 
coincList = []
class event:
    def __init__(self,t,ch,data):
        self.t = t
        self.ch = ch
        self.data = data

class Coincidences:
    def __init__(self):
        self.ch = [0]*16
        self.t = []
        self.data = []
        
    def addEvent(self,ch,t,data):
        self.ch[ch] += 1
        self.time.append(t)
        self.data.append(data)
        
class CoicidenceData:
    def __init__(self,ch,t,data):
        self.ch = [ch]
        self.t = [t]
        self.data = [data]
        
    def addEvent(self,ch,t,data):
        self.ch.append(ch)
        self.t.append(t)
        self.data.append(data)
        
    def writeToDisk(self,filepath,filename):
        
        sortedCh = sorted(self.ch)
        self.t = [x for _,x in sorted(zip(self.ch,self.t),key=lambda pair:pair[0])]
        self.data = [x for _,x in sorted(zip(self.ch,self.data),key=lambda pair:pair[0])]
        self.ch = sortedCh
        
        coincChannelsStr = '_'.join([str(i) for i in sortedCh])
        
        coincfile = f"{filepath}/{filename}_coinc_{coincChannelsStr}.txt" 
        
        if coincfile in openedfiles:
            f = open(coincfile, 'a')
        else:
            coincList.append(sortedCh)
            openedfiles.append(coincfile)
            f = open(coincfile,"w+")
            f.write(f"Data reduction code output.\nEach line saves the information for every channel in the coincidence. Note that currently waveforms are not saved in this data reduction. \nEach channel has the following headers: BOARD;CHANNEL;TIMETAG;ENERGY;ENERGYSHORT;FLAGS;PROBE_CODE;SAMPLES\nNumber of channels in this coincidence: {len(self.ch)}\n")
        
        f.write(';'.join(self.data) + '\n')
        f.close()
        
     
def ReadInData(filepath, filename, savefilepath, coincWindow):
    counter = 0

    with open(f"{filepath}/{filename}.CSV") as f:
        next(f)  # skip header

        for i, line in enumerate(f):
            data = line.strip().split(';')
            
            # Only keep the first 6 fields, ignoring waveform samples beyond that
            cleanLine = ';'.join(data[:6])

            t  = int(data[2]) / 1e3
            ch = int(data[1])

            if counter == 0:
                coinc = CoicidenceData(ch, t, cleanLine)
            else:
                dt = np.absolute(coinc.t[0] - t)

                if dt <= coincWindow:
                    coinc.addEvent(ch, t, cleanLine)
                else:
                    coinc.writeToDisk(filepath=savefilepath, filename=filename)
                    coinc = CoicidenceData(ch, t, cleanLine)

            counter += 1
            if counter % 100000 == 0:
                print(f'{counter} lines sorted')        
                


file = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_06_24/2026_06_24_Daily_LSC_calibration_bck_no_coinc_new_settings_lower_thresh/RAW/SDataR_2026_06_24_Daily_LSC_calibration_bck_no_coinc_new_settings_lower_thresh.CSV')
# filename = 'SDataR_NaI_annulus_LS_Cs137_NaI_1_2_3_4_triple_coinc.CSV'

coincWindow = 0 #ns


filename = str(file).split('.CSV')[0].split('/')[-1]
filepath = str(file).split('.CSV')[0].split(filename)[0]

savefilepath = f"{filepath}/coinc_sorted_{coincWindow}ns"
Path(savefilepath).mkdir(parents=True, exist_ok = True)




ReadInData(filepath=filepath,filename=filename,savefilepath=savefilepath,coincWindow=coincWindow)

print(coincList)

