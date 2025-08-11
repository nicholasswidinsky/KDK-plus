import pandas as pd
from thresh_coinc_fit import erf_gauss_fit
from temperature_run import temperature_reader
import os

directory_str = 'LSC_campaign'
directory = os.fsencode(directory_str)
name_run = []

# put all the file name in a list
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    count = directory_str + "\\" + filename
    if len(os.listdir(count)) > 5:
        name_run.append(filename)
print(name_run)

count = 1
L = []
for file in name_run:
    inpath = '\\Boulot\\Pycharm\\Projects\\KDK+\\LSC_campaign\\' + file + '\\'

    # show where we are
    print(str(count) + "/" + str(len(name_run)))

    if os.path.exists(inpath + 'RAW\\SDataR_' + file + '_coinSorted.csv'):
        dfc = pd.read_csv(inpath + 'RAW\\SDataR_' + file + '_coinSorted.csv', delimiter='\t')
        erf_gauss_fit(dfc, csv=L, plot=False, name=file)
        print("Sample " + file + " done")

        temp_file = "\\Boulot\\Pycharm\\Projects\\KDK+\\Temperature\\LSC_campaign\\temperature_" + file + ".csv"
        if not os.path.exists(temp_file):
            print("No temperature measurement for the sample : " + file)
        else:
            T = temperature_reader(temp_file)
            L[-1]['Temp_mean'] = T[0]
            L[-1]['Temp_strd_dev'] = T[1]
    else:
        print(file + ' is not usable')
    count += 1

data = pd.DataFrame.from_dict(L)
data.to_csv('Figures\\LSC_campaign\\data_analysis\\data.csv', sep=";")
print("Done !!!!")
