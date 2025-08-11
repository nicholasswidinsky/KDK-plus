#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 17:36:58 2024

@author: arnaudlemaire
"""


## Code pour analyser les données des expériences de l'annulus
## Arnaud lemaire 03/05/2024
## Fonctionne avec coinsorted csv files


import os
import numpy as np
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from math import exp
from math import pi
from math import sqrt
from scipy import stats

import pathlib # For paths

import time # For timing
import inspect # To get function names
from deprecated import deprecated # Not sure this works

import subprocess # FOr c++ calls

import collections # For bar charts of text data

from scipy.stats import poisson, norm, truncnorm

import pandas as pd

from uncertainties import ufloat

from num2tex import num2tex

import xml.etree.ElementTree as ET # Parsing xml

from scipy.optimize import curve_fit
 

## Tableau contenant les chemins des données à afficher


def figure(L,title, xmin, xmax,log,nb_bins,d=False):
    n=len(L)
    dfl=[]
    bins=np.linspace(xmin,xmax,nb_bins)
    plt.subplots(figsize=(16, 9))
    nb_events= ''
    for i in range(0,n):
        dfl= pd.read_csv(root+str(L[i][0]), delimiter= '\t' )
        for j in range(0,len(L[i][1])):
            if L[i][1][j][0] == 'energy':
                if len(L[i][1][j][1])==1:
                    dfc = coincidence2(dfl , L[i][1][j][2] )
                    plt.hist(dfc['E'+str(L[i][1][j][1][0])+'_AU'],histtype='step',density=d,label=L[i][1][j][3],color=L[i][1][j][4], bins=bins)
                    plt.xlabel('E'+str(L[i][1][j][1][0])+'_AU')
                    nb_events += str(L[i][1][j][4])+' curve = '+ str(len(dfc)) + ' events \n'
                else:
                    raise NameError('OneChannelOnly')
            if L[i][1][j][0] =='time':
                dfc = coincidence2(dfl , L[i][1][j][2] )
                ps_ns = 1e-3
                print(ps_ns*(dfc['T'+str(L[i][1][j][1][0])+'_ps'] - dfc['T'+str(L[i][1][j][1][1])+'_ps']))
                dft=ps_ns*(dfc['T'+str(L[i][1][j][1][0])+'_ps'] - dfc['T'+str(L[i][1][j][1][1])+'_ps'])
                plt.hist(dft,histtype='step',density=d,label=L[i][1][j][3],color=L[i][1][j][4], bins=bins)
                plt.xlabel('T'+str(L[i][1][j][1][0])+' - T'+str(L[i][1][j][1][1])+'_ns')
                nb_events += str(L[i][1][j][4])+' curve = '+ str(len(dft)) + ' events \n'

    if log==True:
        plt.yscale('log')
    plt.xlim(xmin, xmax)
    paths = ''
    for el in L:
        paths += '& ' + str(el[0]) +'\n'
    plt.title(str(title) + '\n' + root + '\n'+ paths + '\n'+ nb_events +'bins=' + str(len(bins)))

    plt.legend(loc='best',fontsize= 'small'); plt.ylabel('Counts')
#    plt.savefig(root + str(title),dpi='figure')
    return plt.show()


    
## Fonction coincidence qui à partir d'un dataframe de coinsorted.csv, renvoie un 
## dataframe composé des éléments coincidents de certains canaux 

def coincidence(df,Ch_coincidence):
    dfCo=len(Ch_coincidence)*[0]
    k=1
    if len(Ch_coincidence) < 2:
        return df
    else:
        dfCo[0]= df[ df['E'+str(Ch_coincidence[0])+'_AU']>0 ]
        for i in Ch_coincidence[1:]:
            dfCo[k]= dfCo[k-1][dfCo[k-1]['E'+str(i)+'_AU']>0]
            k+=1
        dfCoin = dfCo[-1]
        return dfCoin

## Coincidence 2 intègre un cut en énergie dans le canal coincident (cut en temps à ajouter)

def coincidence2(df,Ch_coincidence):
    dfCo=len(Ch_coincidence)*[0]
    k=1
    if len(Ch_coincidence) < 2:
        return df
    else:
        dfCo[0]= df[ df['E'+str(Ch_coincidence[0][0])+'_AU']>0 ]
        for el in Ch_coincidence[1:]:
            if len(el) == 1:
                dfCo[k]= dfCo[k-1][dfCo[k-1]['E'+str(el[0])+'_AU']>0]
                k+=1
            else:
                condition_energy = ( dfCo[k-1]['E'+str(el[0])+'_AU'] > el[1][0]) & (dfCo[k-1]['E'+str(el[0])+'_AU'] < el[1][1])
                dfCo[k]= dfCo[k-1].loc[condition_energy,:].copy()
                k+=1
        dfCoin = dfCo[-1]
        return dfCoin

def gauss(x,a,b,m,sigma):
    return (a/(sigma*sqrt(2*pi)))*np.exp(-((x-m)**2)/(2*(sigma**2))) + b 


## Fonction pour ajouter un fit gaussien sur un histogram

# Utilisation de la fonction norm.pdf et stats.norm.fit pour ajuster une gaussiènne sur laleau lui même
# fonction normlisée donc la hauteur n'est pas ajustée avec l'histogramme

def fit(F,title, xmin, xmax,log,nb_bins,d=True):
    n=len(F)
    dfl=[]
    bins=np.linspace(xmin,xmax,nb_bins)
    plt.subplots(figsize=(16, 9))
    couleur=['y--','m--','c--','b--','g--', 'r--']
    nb_events= ''
    for i in range(0,n):
        dfl= pd.read_csv(root+str(F[i][0]), delimiter= '\t' )
        for j in range(0,len(F[i][1])):
            if F[i][1][j][0] == 'energy':
                dfc = coincidence2(dfl , F[i][1][j][2] )

                condition_energy = ( dfc['E'+str(F[i][1][j][1])+'_AU'] > F[i][1][j][3][0]) & (dfc['E'+str(F[i][1][j][1])+'_AU'] < F[i][1][j][3][1])
                df_fit= dfc.loc[condition_energy,:].copy()
                params = stats.norm.fit(df_fit['E'+str(F[i][1][j][1])+'_AU'])
                mean = params[0]
                fwhm = 2*np.sqrt(2*np.log(2))*params[1]
                
                binsfit=np.linspace(F[i][1][j][3][0],F[i][1][j][3][1],100)
                f = (F[i][1][j][3][1] - F[i][1][j][3][0]) / (xmax - xmin)
                plt.plot(binsfit, f*norm.pdf(binsfit,params[0],params[1]), couleur.pop(),label='fit: m= ' + str(round(mean)) + ' fwhm= ' + str(round(fwhm)))
                
                plt.hist(dfc['E'+str(F[i][1][j][1])+'_AU'],histtype='step',density=d,label=F[i][1][j][4],color=F[i][1][j][5], bins=nb_bins)
                plt.xlabel('E'+str(F[i][1][j][1])+'_AU')
                nb_events += str(F[i][1][j][5])+' curve = '+ str(len(dfc)) + ' events \n'

            if F[i][1][j][0] =='time':
                dfc = coincidence2(dfl , F[i][1][j][2] )
                ps_ns = 1e-3

                dft=ps_ns*(dfc['T'+str(F[i][1][j][1][0])+'_ps'] - dfc['T'+str(F[i][1][j][1][1])+'_ps'])
                print(dft)
                plt.hist(dft,histtype='step',density=d,label=F[i][1][j][4],color=F[i][1][j][5], bins=bins)
                plt.xlabel('T'+str(F[i][1][j][1][0])+' - T'+str(F[i][1][j][1][1])+'_ns')    
                nb_events += str(F[i][1][j][5])+' curve = '+ str(len(dft)) + ' events \n'
                
                condition_time = (dft > F[i][1][j][3][0]) & (dft < F[i][1][j][3][1])
                df_fit= dft.loc[condition_time].copy()
                params = stats.norm.fit(df_fit)
                mean = params[0]
                fwhm = 2*np.sqrt(2*np.log(2))*params[1]
                sigma= params[1]
                print(sigma)
                print(df_fit)
                binsfit=np.linspace(F[i][1][j][3][0],F[i][1][j][3][1],60)

                plt.plot(binsfit,norm.pdf(binsfit,params[0],params[1]), couleur.pop(),label='fit: m= ' + str(round(mean)) + ' fwhm= ' + str(round(fwhm)))
                
    if log==True:
        plt.yscale('log')
    plt.xlim(xmin, xmax)
    paths = ''
    for el in F:
        paths += '& ' + str(el[0]) +'\n'
    plt.title(str(title) + '\n' + root + '\n'+ paths + '\n'+ nb_events +'bins=' + str(len(bins)))
    plt.legend(loc='best',fontsize= 'medium'); plt.ylabel('Counts_density')
#    plt.savefig(root + str(title),dpi='figure')
    return plt.show()

    
#%%

if __name__ == '__main__':

# Liste L de la forme: L =[ [path,channels] ,  ...]
# Avec channels= [[energy/time/flag,channel,coincidence,title,color], ...]
    
    
    root = '/Volumes/KDK+_Arnaud/KDK+/teflon_vessel/annulus2/'
    
    if True:    
        path1='le_source_above/RAW/SDataR_le_source_above_coinSorted.csv'
        path2='cfd_source_above_param_100p_150ns/RAW/SDataR_cfd_source_above_param_100p_150ns_coinSorted.csv'
        path3='cfd_source_above_param_100p_200ns/RAW/SDataR_cfd_source_above_param_100p_200ns_coinSorted.csv'
        path4='cfd_source_above_param_100p_300ns/RAW/SDataR_cfd_source_above_param_100p_300ns_coinSorted.csv'
        path5='cfd_source_above_param_100p_250ns/RAW/SDataR_cfd_source_above_param_100p_250ns_coinSorted.csv'
        channels1=[['energy',[4],[0,2,4] ,'Leading_Edge_coincident','blue'],
                   ['energy',[4],[] ,'Leading_Edge_total','darkblue']]
        channels2=[['energy',[4],[0,2,4] ,'CFD_100p_150ns_coincident','red'],
                   ['energy',[4],[] ,'CFD_100p_150ns_total','darkred']]
        channels3=[['energy',[4],[0,2,4] ,'CFD_100p_200ns_coincident','green'],
                   ['energy',[4],[] ,'CFD_100p_200ns_total','darkgreen']]
        channels4=[['energy',[4],[0,2,4] ,'CFD_100p_300ns_coincident','magenta'],
                   ['energy',[4],[] ,'CFD_100p_300ns_total','darkmagenta']]
        channels5=[['energy',[4],[0,2,4] ,'CFD_100p_250ns_coincident','peru'],
                   ['energy',[4],[] ,'CFD_100p_250ns_total','saddlebrown']]
        L=[[path1,channels1],[path2,channels2],[path3,channels3],[path4,channels4],[path5,channels5]]
    if False:    
        path1='le_source_above/RAW/SDataR_le_source_above_coinSorted.csv'
        path2='cfd_source_above_param_100p_150ns/RAW/SDataR_cfd_source_above_param_100p_150ns_coinSorted.csv'
        path3='cfd_source_above_param_100p_200ns/RAW/SDataR_cfd_source_above_param_100p_200ns_coinSorted.csv'
        path4='cfd_source_above_param_100p_300ns/RAW/SDataR_cfd_source_above_param_100p_300ns_coinSorted.csv'
        path5='cfd_source_above_param_100p_250ns/RAW/SDataR_cfd_source_above_param_100p_250ns_coinSorted.csv'
        channels1=[
                   ['time',[0,4],[0,2,4] ,'Leading_Edge_total','darkblue']]
        channels2=[
                   ['time',[0,4],[0,2,4] ,'CFD_100p_150ns_total','darkred']]
        channels3=[
                   ['time',[0,4],[0,2,4] ,'CFD_100p_200ns_total','darkgreen']]
        channels4=[
                   ['time',[0,4],[0,2,4] ,'CFD_100p_300ns_total','darkmagenta']]
        channels5=[
                   ['time',[0,4],[0,2,4] ,'CFD_100p_250ns_total','saddlebrown']]
        L=[[path1,channels1],[path2,channels2],[path3,channels3],[path4,channels4],[path5,channels5]]
    
    if False:
        L=[['le_source_above/RAW/SDataR_le_source_above_coinSorted.csv',[['time',[0,2],[0,2,4] ,'LSC_coincidence_Leading_Edge','black'],
        ['time',[0,4],[0,2,4] ,'LSC_annulus_coincidence_Leading_Edge','blue']]],
       ['cfd_source_above/RAW/SDataR_cfd_source_above_coinSorted.csv',[['time',[0,2],[0,2,4] ,'LSC_coincidence_CFD','red'],
        ['time',[0,4],[0,2,4] ,'LSC_annulus_coincident_CFD','orange']]]
       ]

    if True: figure(L,'NaI_spectra_annulus_coincidence_Teflon_vessel_65Zn_above_triggering_parameter' ,0,2500,True,200)
    
    plt.clf()
    
#%%

    root = '/Volumes/KDK+_Arnaud/KDK+/teflon_vessel/annulus2/'
# Liste F de la forme: F =[ [path,channels] ,  ...]
# Avec channels= [[energy/time/flag,channel,coincidence,fit interval, title,color], ...]
    if True:    
        path1='le_source_above/RAW/SDataR_le_source_above_coinSorted.csv'
        path2='cfd_source_above_param_100p_150ns/RAW/SDataR_cfd_source_above_param_100p_150ns_coinSorted.csv'
        path3='cfd_source_above_param_100p_200ns/RAW/SDataR_cfd_source_above_param_100p_200ns_coinSorted.csv'
        path4='cfd_source_above_param_100p_300ns/RAW/SDataR_cfd_source_above_param_100p_300ns_coinSorted.csv'
        path5='cfd_source_above_param_100p_250ns/RAW/SDataR_cfd_source_above_param_100p_250ns_coinSorted.csv'
        channels1=[
                   ['time',[0,4],[0,2,4],[-100,-60] ,'Leading_Edge_total','darkblue']]
        channels2=[
                   ['time',[0,4],[0,2,4],[-430,-350] ,'CFD_100p_150ns_total','darkred']]
        channels3=[
                   ['time',[0,4],[0,2,4],[-490,-380] ,'CFD_100p_200ns_total','darkgreen']]
        channels4=[
                   ['time',[0,4],[0,2,4],[-550,-450] ,'CFD_100p_300ns_total','darkmagenta']]
        channels5=[
                   ['time',[0,4],[0,2,4],[-495,-425] ,'CFD_100p_250ns_total','saddlebrown']]
        F=[[path1,channels1],[path2,channels2],[path3,channels3],[path4,channels4],[path5,channels5]]
        
    if False:    
        path1='le_source_above/RAW/SDataR_le_source_above_coinSorted.csv'
        path2='cfd_source_above/RAW/SDataR_cfd_source_above_coinSorted.csv'
#        path3=
#        path4=
        channels1=[['time',[0,2],[0,2,4],[-8,2],'LSC_coincidence_Leading_Edge','black']]
        channels2=[['time',[0,2],[0,2,4],[-8,2] ,'LSC_coincidence_CFD','red'],
                   ['time',[0,4],[0,2,4],[-375,-335] ,'LSC_annulus_coincident_CFD','orange']]
        F=[[path2,channels2]]
        
    if True: fit(F,'CoincidenceTime_NaI_annulus_Teflon_vessel_65Zn_above_triggering_parameter' ,-600,0,False,200)
    
    plt.clf()
    
    #%%
    #Test du crystal BICRON

    root = '/Volumes/KDK+_Arnaud/KDK+/time_resolution_test/DAQ/'
# Liste F de la forme: F =[ [path,channels] ,  ...]
# Avec channels= [[energy/time/flag,channel,coincidence,fit interval, title,color], ...]
    if True:    
        path1='bicronNaI_65Zn_teflon_coincidence_LE/RAW/SDataR_bicronNaI_65Zn_teflon_coincidence_LE_coinSorted.csv'
#        path2='cfd_source_above_param_100p_150ns/RAW/SDataR_cfd_source_above_param_100p_150ns_coinSorted.csv'
#        path3='cfd_source_above_param_100p_200ns/RAW/SDataR_cfd_source_above_param_100p_200ns_coinSorted.csv'

        channels1=[['time',[0,4],[0,2,4],[-230,-200] ,'LSC-NaI','black'],
                   ['time',[0,2],[0,2,4],[-15,5] ,'LSC-LSC','red']]
#        channels2=[
#                   ['time',[0,4],[0,2,4],[-430,-350] ,'CFD_100p_150ns_total','darkred']]
#        channels3=[
#                   ['time',[0,4],[0,2,4],[-490,-380] ,'CFD_100p_200ns_total','darkgreen']]

        F=[[path1,channels1]]
        
    if False:    
        path1='bicronNaI_65Zn_teflon_coincidence_LE/RAW/SDataR_bicronNaI_65Zn_teflon_coincidence_LE_coinSorted.csv'
#        path2='cfd_source_above/RAW/SDataR_cfd_source_above_coinSorted.csv'
#        path3=
#        path4=
        channels1=[['time',[0,4],[0,2,4],'LSC-NaI','black'],
                   ['time',[0,2],[0,2,4],'LSC-LSC','red']]
#        channels2=[['time',[0,2],[0,2,4],[-8,2] ,'LSC_coincidence_CFD','red'],
#                   ['time',[0,4],[0,2,4],[-375,-335] ,'LSC_annulus_coincident_CFD','orange']]
        L=[[path1,channels1]]
        
    if True: fit(F,'CoincidenceTime_NaI_annulus_Teflon_vessel_65Zn_above_triggering_parameter' ,-250,30,False,200)

    if False: figure(L,'Bicron_NaI_coincidence_Teflon_vessel_65Zn_Time_resolution' ,-300,20,True,200)

    plt.clf()
    
    
        #%%
        ## Code pour la comparaison des LSC

    root = '/Volumes/KDK+_Arnaud/KDK+/LSC_campaign/DAQ/'
# Liste F de la forme: F =[ [path,channels] ,  ...]
# Avec channels= [[energy/time/flag,channel,coincidence ,fit interval, title,color], ...]
    if True:    
        path1='UG_pure/RAW/SDataR_UG_pure_coinSorted.csv'
        path2='LLT_pure/RAW/SDataR_LLT_pure_coinSorted.csv'
        path3='LLT_3ml_1M/RAW/SDataR_LLT_3ml_1M_coinSorted.csv'

        # channels1=[['time',[0,4],[0,2,4],[-230,-200] ,'LSC-NaI','black'],
        #            ['time',[0,2],[0,2,4],[-15,5] ,'LSC-LSC','red']]

        # channels1=[['energy',[0],[[0],[2,[1700,2700]]] ,'60°','red'],
        #            ['energy',[0],[[0],[2]], '60°_nocut','darkred']]
        # channels1=[['energy',[0],[[0],[4,[1300,2000]]] ,'90°','blue'],
        #            ['energy',[0],[[0],[4]], '90°_nocut','darkblue']]
        channels1=[['energy',[0],[[0],[2]] ,'60°_UG','blue'],
                    ['energy',[0],[[0],[4]], '90°_UG','red']]
        channels2=[['energy',[0],[[0],[2]] ,'60°_LLT','darkblue'],
                    ['energy',[0],[[0],[4]], '90°_LLT','darkred']]
        channels3=[['energy',[0],[[0],[2]] ,'60°_LLT_1M_3mL','blue'],
                    ['energy',[0],[[0],[4]], '90°_LLT_1M_3mL','red']]
        L=[[path3,channels3],[path2,channels2]]

# Liste L de la forme: L =[ [path,channels] ,  ...]
# Avec channels= [[energy/time/flag,channel,coincidence,title,color], ...]
# avec coincidence= [[channel,cut],..]

    if False:    
        path1='UG_pure/RAW/SDataR_UG_pure_coinSorted.csv'
        path2='LLT_pure/RAW/SDataR_LLT_pure_coinSorted.csv'

#        channels1=[['energy',[0],[0,2],'60°','black'],
#                   ['energy',[0],[0,4],'90°','red']]
        channels1=[['energy',[0],[[0],[2,[1700,2700]]],[1100,1900] ,'60°','red'],
                   ['energy',[0],[[0],[4,[1400,1900]]],[1600,2400] ,'90°','orange'],
                   ]
        F=[[path1,channels1]]
        
    if False: fit(F,'Calibration_UG_pure' ,0,3000,False,100)

    if True: figure(L,'Calibration_UG_pure' ,0,3200,False,100)

    plt.clf()
    