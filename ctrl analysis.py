# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:49:10 2021

@author: wiesbrock
"""

#Control analysis

import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
import os
import scipy.stats as stats
import seaborn as sns
import scikit_posthocs as sp

#Zeitslots

def crosscorr(datax, datay, lag=0):
    """ Lag-N cross correlation. 
    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length

    Returns
    ----------
    crosscorr : float
    """
    return datax.corr(datay.shift(lag))


folder=r'C:\Users\wiesbrock\Desktop\Analysen\Ctrl'
folder_list=folder+'\\*'
folder_list=glob.glob(folder_list)

os.chdir(folder)
try:
    os.mkdir('plots')
except:
    print('Folder already exists')
    
peak_list=[]
std_list=[]
dist_peaks=[]

# window 1
for i in folder_list:
    files=str(i+'\\video\\results\\manual\\time_series.xlsx')
    df =pd.read_excel(files,sheet_name=None)
    a=df['detrended']
    header=list(a.columns.values)
    names=header
    corr_matrix=np.zeros((len(names),len(names)))
    index_max_corr=np.zeros((len(names),len(names)))
    b=np.zeros((40,1))
    peaks_all=np.zeros((len(names),1))
    print(len(names))
    rand=np.random.randint(0,(len(a)-100),1)
    rand=rand.astype(int)
    start=rand[0]
    stop=start+100
    c=a[start:stop]
    
    for j in range(len(names)):
        
        zscore=stats.zscore(c[names[j]])
        peaks=np.where(zscore>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks=len(peak_diff[peak_diff>=10])
        peak_list.append(number_of_peaks)
        plt.figure(1)
        plt.xlim(-3,3)
        hist,edges=np.histogram(zscore)
        hist=hist/np.max(hist)
        plt.ylabel('Amount normalized to max')
        plt.xlabel('Z score')
        plt.plot(edges[:-1],hist)
        dist_peaks.append(edges[np.where(hist==np.max(hist))][0])
   
        for k in range(len(names)):
            for n in range(40):
                trace_1=c[names[j]]
                trace_2=c[names[k]]
                b[n]=crosscorr(trace_1,trace_2, lag=n)
                corr_matrix[j,k]=np.max(b)
                index_max_corr[j,k]=np.where(b==np.max(b))[0]
                
    win1=corr_matrix
                
                #Crosscorr für verschiedene Windows
'''            
plt.figure()       
sns.boxplot(data=dist_peaks, width=1)
plt.plot(np.linspace(-1,1,3),np.linspace(1.96,1.96,3),'k--')
plt.plot(np.linspace(-1,1,3),np.linspace(-1.96,-1.96,3),'k--')
plt.xticks([])
plt.ylabel('Peak histogram z-score')
plt.xlabel('Dashed lines indicate Significance boundary for alpha=0.05')
'''

# window 2
for i in folder_list:
    files=str(i+'\\video\\results\\manual\\time_series.xlsx')
    df =pd.read_excel(files,sheet_name=None)
    a=df['detrended']
    header=list(a.columns.values)
    names=header
    corr_matrix=np.zeros((len(names),len(names)))
    index_max_corr=np.zeros((len(names),len(names)))
    b=np.zeros((40,1))
    peaks_all=np.zeros((len(names),1))
    rand=np.random.randint(0,(len(a)-100),1)
    rand=rand.astype(int)
    start=rand[0]
    stop=start+100
    c=a[start:stop]
    
    for j in range(len(names)):
        
        zscore=stats.zscore(c[names[j]])
        peaks=np.where(zscore>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks=len(peak_diff[peak_diff>=10])
        peak_list.append(number_of_peaks)
        plt.figure(1)
        plt.xlim(-3,3)
        hist,edges=np.histogram(zscore)
        hist=hist/np.max(hist)
        plt.ylabel('Amount normalized to max')
        plt.xlabel('Z score')
        plt.plot(edges[:-1],hist)
        dist_peaks.append(edges[np.where(hist==np.max(hist))][0])
   
        for k in range(len(names)):
            for n in range(40):
                trace_1=c[names[j]]
                trace_2=c[names[k]]
                b[n]=crosscorr(trace_1,trace_2, lag=n)
                corr_matrix[j,k]=np.max(b)
                index_max_corr[j,k]=np.where(b==np.max(b))[0]
                
        win2=corr_matrix
                
                #Crosscorr für verschiedene Windows
# window 3
for i in folder_list:
    files=str(i+'\\video\\results\\manual\\time_series.xlsx')
    df =pd.read_excel(files,sheet_name=None)
    a=df['detrended']
    header=list(a.columns.values)
    names=header
    corr_matrix=np.zeros((len(names),len(names)))
    index_max_corr=np.zeros((len(names),len(names)))
    b=np.zeros((40,1))
    peaks_all=np.zeros((len(names),1))
    rand=np.random.randint(0,(len(a)-100),1)
    rand=rand.astype(int)
    start=rand[0]
    stop=start+100
    c=a[start:stop]
    
    for j in range(len(names)):
        
        zscore=stats.zscore(c[names[j]])
        peaks=np.where(zscore>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks=len(peak_diff[peak_diff>=10])
        peak_list.append(number_of_peaks)
        plt.figure(1)
        plt.xlim(-3,3)
        hist,edges=np.histogram(zscore)
        hist=hist/np.max(hist)
        plt.ylabel('Amount normalized to max')
        plt.xlabel('Z score')
        plt.plot(edges[:-1],hist)
        dist_peaks.append(edges[np.where(hist==np.max(hist))][0])
   
        for k in range(len(names)):
            for n in range(40):
                trace_1=c[names[j]]
                trace_2=c[names[k]]
                b[n]=crosscorr(trace_1,trace_2, lag=n)
                corr_matrix[j,k]=np.max(b)
                index_max_corr[j,k]=np.where(b==np.max(b))[0]
                
        win3=corr_matrix
                
# window 4
for i in folder_list:
    files=str(i+'\\video\\results\\manual\\time_series.xlsx')
    df =pd.read_excel(files,sheet_name=None)
    a=df['detrended']
    header=list(a.columns.values)
    names=header
    corr_matrix=np.zeros((len(names),len(names)))
    index_max_corr=np.zeros((len(names),len(names)))
    b=np.zeros((40,1))
    peaks_all=np.zeros((len(names),1))
    rand=np.random.randint(0,(len(a)-100),1)
    rand=rand.astype(int)
    start=rand[0]
    stop=start+100
    c=a[start:stop]
    
    for j in range(len(names)):
        
        zscore=stats.zscore(c[names[j]])
        peaks=np.where(zscore>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks=len(peak_diff[peak_diff>=10])
        peak_list.append(number_of_peaks)
        plt.figure(1)
        plt.xlim(-3,3)
        hist,edges=np.histogram(zscore)
        hist=hist/np.max(hist)
        plt.ylabel('Amount normalized to max')
        plt.xlabel('Z score')
        plt.plot(edges[:-1],hist)
        dist_peaks.append(edges[np.where(hist==np.max(hist))][0])
   
        for k in range(len(names)):
            for n in range(40):
                trace_1=c[names[j]]
                trace_2=c[names[k]]
                b[n]=crosscorr(trace_1,trace_2, lag=n)
                corr_matrix[j,k]=np.max(b)
                index_max_corr[j,k]=np.where(b==np.max(b))[0]
                
    win4=corr_matrix
                
# window 5
for i in folder_list:
    files=str(i+'\\video\\results\\manual\\time_series.xlsx')
    df =pd.read_excel(files,sheet_name=None)
    a=df['detrended']
    header=list(a.columns.values)
    names=header
    corr_matrix=np.zeros((len(names),len(names)))
    index_max_corr=np.zeros((len(names),len(names)))
    b=np.zeros((40,1))
    peaks_all=np.zeros((len(names),1))
    print(len(names))
    rand=np.random.randint(0,(len(a)-100),1)
    rand=rand.astype(int)
    start=rand[0]
    stop=start+100
    c=a[start:stop]
    
    for j in range(len(names)):
        
        zscore=stats.zscore(c[names[j]])
        peaks=np.where(zscore>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks=len(peak_diff[peak_diff>=10])
        peak_list.append(number_of_peaks)
        plt.figure(1)
        plt.xlim(-3,3)
        hist,edges=np.histogram(zscore)
        hist=hist/np.max(hist)
        plt.ylabel('Amount normalized to max')
        plt.xlabel('Z score')
        plt.plot(edges[:-1],hist)
        dist_peaks.append(edges[np.where(hist==np.max(hist))][0])
   
        for k in range(len(names)):
            for n in range(40):
                trace_1=c[names[j]]
                trace_2=c[names[k]]
                b[n]=crosscorr(trace_1,trace_2, lag=n)
                corr_matrix[j,k]=np.max(b)
                index_max_corr[j,k]=np.where(b==np.max(b))[0]
                
    win5=corr_matrix
    
import seaborn as sns
sns.set_palette("rocket_r")
plt.figure()
plot_data=win1,win2,win3,win4,win5
sns.violinplot(data=plot_data).set_title('Crosscorrelation')
#plt.sca(f[ax_swarm])
#plt.xticks((0,1,2,3,4,5), ['Pre1','Pre2','Stim','Post1','Post2', 'Ctrl'])
plt.xticks(range(5),['Win1','Win2','Win3','Win4','Win5'])

                
                

