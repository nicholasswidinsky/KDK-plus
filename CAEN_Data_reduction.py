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
        
        coincChannelsStr = '_'.join([str(i) for i in sortedCh])
        
        coincfile = f"{filepath}/{filename}_coinc_{coincChannelsStr}.txt" 
        
        if coincfile in openedfiles:
            f = open(coincfile, 'a')
        else:
            coincList.append(sortedCh)
            openedfiles.append(coincfile)
            f = open(coincfile,"w+")
            f.write(f"Data reduction code output.\nEach line saves the information for every channel in the coincidence. Note that currently waveforms are not saved in this data reduction. \nEach channel has the following headers: BOARD;CHANNEL;TIMETAG;ENERGY;ENERGYSHORT;FLAGS;PROBE_CODE;SAMPLES\nNumber of channels in this coincidenc: {len(self.ch)}\n")
        
        f.write(';'.join(self.data) + '\n')
        f.close()
        
        

def ReadInData(filepath,filename,savefilepath,coincWindow,):
    # coinc_list = []
    # current_coinc = Coincidences()
    # first_event_time = None
    counter = 0

    with open(f"{filepath}/{filename}.CSV") as f:
        next(f)

        for i, line in enumerate(f):
            
            lines = line.split('\n')[0]
            # print(lines)
            if counter == 0:
               
                data = line.strip().split(';')

                t = int(data[2]) / 1e3
                ch = int(data[1])

                coinc = CoicidenceData(ch,t,lines)
            
            else: 
                data = line.strip().split(';')

                t = int(data[2]) / 1e3
                ch = int(data[1])
                
                dt = np.absolute(coinc.t[0] - t)
                
                if dt <= coincWindow:
                    coinc.addEvent(ch,t,lines)
                    
                else:
                    coinc.writeToDisk(filepath=savefilepath,filename=filename)
                    
                    coinc = CoicidenceData(ch,t,lines)
            counter += 1 
            if counter % 100000 == 0:
                print(f'{counter} lines sorted')
            # counter += 1
            # # if counter % 100000 == 0:
            # #     print(f"Read in {counter} events")

            # if first_event_time is None:
            #     first_event_time = t

            # dt = abs(first_event_time - t)

            # if dt > CoincWindow:
            #     if current_coinc.Channels not in CoincChannels:
            #         del current_coinc
            #     else:
            #         coinc_list.append(current_coinc)
            #     current_coinc = Coincidence()
            #     first_event_time = t

            # current_coinc.AddEvent(evt, ch)

        # coinc_list.append(current_coinc)


file = '/home/nick/PhD/KDK+/Annulus_Compton_scatter_V1/2026_01_28/NaI_annulus_LS_2_Cs137_All_NaI_higher_HV/RAW/SDataR_NaI_annulus_LS_2_Cs137_All_NaI_higher_HV.CSV'
# filename = 'SDataR_NaI_annulus_LS_Cs137_NaI_1_2_3_4_triple_coinc.CSV'

filename = file.split('.CSV')[0].split('/')[-1]
filepath = file.split('.CSV')[0].split(filename)[0]

savefilepath = f"{filepath}/coinc_sorted"
Path(savefilepath).mkdir(parents=True, exist_ok = True)

coincWindow = 500 #ns

ReadInData(filepath=filepath,filename=filename,savefilepath=savefilepath,coincWindow=coincWindow)

print(coincList)

