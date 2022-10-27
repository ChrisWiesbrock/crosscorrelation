# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 15:01:59 2022
@author: wiesbrock
"""

from IPython import get_ipython
#get_ipython().magic('reset -sf')

import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
import os
import scipy.stats as stats
import seaborn as sns
import scikit_posthocs as sp




folder=r'C:\Users\wiesbrock\Desktop\Daten_corr_analyse\190825_5.2\video'
search=folder+'*'
folder_list=glob.glob(search)
os.chdir(folder)
try:
    os.mkdir('plots')
    os.mkdir('plots\\single traces')
    os.mkdir('plots\\Control_single traces')
    os.mkdir('plots\\auto')
except:
    print('Folder already exists')

control=r'C:\Users\wiesbrock\Desktop\Daten_corr_analyse\ctrl\181126_5.2\video'
control_search=folder+'*'
control_folder_list=glob.glob(search)




def moving_average(data_set, periods=3):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, mode='valid')

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


for m in range(len(folder_list)):

    #files=str(folder_list[m])+'\\results\\manual\\time_series.xlsx'
    #meta=str(folder_list[m])+'\\Intervals.txt'
    
    files=str(folder+'\\results\\manual\\time_series.xlsx')
    meta=str(folder+'\\Intervals.txt')



    df =pd.read_excel(files,sheet_name=None)

    intervall = pd.read_csv(meta, header = None)
    if np.shape(intervall)==(4,4):
        if int(intervall[2][2])==1:
            intervall[2][2]=int(intervall[2][2])+int(intervall[3][1])
            intervall[3][2]=int(intervall[3][2])+int(intervall[3][1])
            intervall[2][3]=int(intervall[2][3])+int(intervall[3][2])
            intervall[3][3]=int(intervall[3][3])+int(intervall[3][2])
    a=df['detrended']
    stimulus_length=int(intervall[3][2])-int(intervall[2][2])-1
    header=list(a.columns.values)
    activity_overview=np.zeros((5,np.shape(a)[1]))
    v=0
    peaks_all=np.zeros((len(header),1))
    j=0
    for i in header:
        
        zscore=stats.zscore(a[i])
        peaks=np.where(zscore>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks=len(peak_diff[peak_diff>=15])
        peaks_all[v]=number_of_peaks+1
        plt.figure()
        plt.title(files[30:]+' '+i+ 'peaks='+str(number_of_peaks+1))
        #plt.axis('off')
        ax=plt.subplot()
        ax.plot(a[i], 'k')
        ax.plot(a[i][peaks],'r.')
        if np.shape(intervall)==(4,4):
            inter=np.linspace(int(intervall[2][2]),int(intervall[3][2]),len(a[i]))
            max_trace=np.linspace(np.max(a[i])+0.05,np.max(a[i])+0.05,len(a[i]))
            ax.plot(inter,max_trace,'k')
        plt.savefig(folder+'\\plots\\single traces\\' +i+'.svg')
        v=v+1
        thr=0
        
        
        #pre1 activity
        start=int(intervall[2][1])
        zscore_pre=zscore[0:100]
        peaks=np.where(zscore_pre>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks_pre=len(peak_diff[peak_diff>=10])
        if len(peaks)>0:
            number_of_peaks_pre=number_of_peaks_pre+1
        #print(number_of_peaks_pre)
        if number_of_peaks_pre>thr:
            activity_overview[0,j]=1
        
        
        #pre2 activity
        start=int(intervall[2][1])
        zscore_pre=zscore[100:200]
        peaks=np.where(zscore_pre>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks_pre=len(peak_diff[peak_diff>=10])
        if len(peaks)>0:
            number_of_peaks_pre=number_of_peaks_pre+1
        print(number_of_peaks_pre)
        if number_of_peaks_pre>thr:
            activity_overview[1,j]=1
            
        #stim activity
        start=int(intervall[2][2])
        stop=int(intervall[3][2])
        print(start)
        print(stop)
        zscore_pre=zscore[start:stop]
        peaks=np.where(zscore_pre>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        check=peak_diff
        peak_diff=peak_diff
        number_of_peaks_pre=len(peak_diff[peak_diff>=10])
        if len(peaks)>0:
            number_of_peaks_pre=number_of_peaks_pre+1
        #print(number_of_peaks_pre)
        if number_of_peaks_pre>thr:
            activity_overview[2,j]=1
            
        #poststim1 activity
        start=int(intervall[3][2])
        zscore_pre=zscore[start:start+100]
        peaks=np.where(zscore_pre>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks_pre=len(peak_diff[peak_diff>=10])
        if len(peaks)>0:
            number_of_peaks_pre=number_of_peaks_pre+1
        #print(number_of_peaks_pre)
        if number_of_peaks_pre>thr:
            activity_overview[3,j]=1
            
        #poststim2 activity
        
        start=int(intervall[3][2])
        zscore_pre=zscore[-100:]
        peaks=np.where(zscore_pre>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks_pre=len(peak_diff[peak_diff>=10])
        if len(peaks)>0:
            number_of_peaks_pre=number_of_peaks_pre+1
        print(number_of_peaks_pre)
        if number_of_peaks_pre>thr:
            activity_overview[4,j]=1
            
        j=j+1

active_indices=np.where(peaks_all>=-1)[0]
d=pd.DataFrame()


#prestim1


for i in range(len(active_indices)):
    d[header[active_indices[i]]]=a[header[active_indices[i]]]


    

start=int(intervall[3][1])
pre=d[0:start]
names=list(d.columns)
names=np.array(names)[activity_overview[0,:]==1]
corr_matrx=np.zeros((len(names),len(names)))
index_max_corr=np.zeros((len(names),len(names)))
b=np.zeros((40,1))


for i in range(len(names)):
    for k in range(len(names)):
        for n in range(40):
            trace_1=(pre[names[i]])[-stimulus_length:]
            trace_2=(pre[names[k]])[-stimulus_length:]
            b[n]=crosscorr(trace_1,trace_2, lag=n)
        corr_matrx[i,k]=np.max(b)
        index_max_corr[i,k]=np.where(b==np.max(b))[0]
        
pre_matrix_1=np.reshape(corr_matrx,(len(corr_matrx)**2,1))
pre_matrix_1=pre_matrix_1[pre_matrix_1<1]
#print(index_max_corr[activity_overview_after_check[0,:]==1])

if len(names)>=2:
    plt.figure()
    plt.title('Prestim1')
    sns.heatmap(corr_matrx, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')
    #x=np.linspace(0,len(names),len(names))


    plt.savefig(folder+'\\plots\\prestim_1_corr.svg')

    plt.figure()
    plt.title('Prestim1 Lag')
    sns.heatmap(index_max_corr, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\prestim_1_corr_lag.svg')


#prestim2

for i in range(len(active_indices)):
    d[header[active_indices[i]]]=a[header[active_indices[i]]]
   
start=int(intervall[3][1])
pre=d[100:start+100]
names=list(d.columns)
names=np.array(names)[activity_overview[1,:]==1]
corr_matrx=np.zeros((len(names),len(names)))
index_max_corr=np.zeros((len(names),len(names)))
b=np.zeros((40,1))
for i in range(len(names)):
    for k in range(len(names)):
        for n in range(40):
            trace_1=(pre[names[i]])[-stimulus_length:]
            trace_2=(pre[names[k]])[-stimulus_length:]
            b[n]=crosscorr(trace_1,trace_2, lag=n)
        corr_matrx[i,k]=np.max(b)
        index_max_corr[i,k]=np.where(b==np.max(b))[0]
        
pre_matrix_2=np.reshape(corr_matrx,(len(corr_matrx)**2,1))
pre_matrix_2=pre_matrix_2[pre_matrix_2<1]

if len(names)>=2:
    plt.figure()
    plt.title('Prestim_2')
    sns.heatmap(corr_matrx, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    x=np.linspace(0,len(names)-1,len(names))
    plt.yticks(x,names)
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\prestim_2_corr.svg')

    plt.figure()
    plt.title('Prestim_2 Lag')
    sns.heatmap(index_max_corr, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')
    
plt.savefig(folder+'\\plots\\prestim_2_corr_lag.svg')


#stim

for i in range(len(active_indices)):
    d[header[active_indices[i]]]=a[header[active_indices[i]]]

start=int(intervall[2][2])
stop=int(intervall[2][3])
stim=d[start:stop]
names=list(d.columns)
names=np.array(names)[activity_overview[2,:]==1]
corr_matrx=np.zeros((len(names),len(names)))
index_max_corr=np.zeros((len(names),len(names)))
b=np.zeros((40,1))
for i in range(len(names)):
    for k in range(len(names)):
        for n in range(40):
            b[n]=crosscorr(stim[names[i]],stim[names[k]], lag=n)
        corr_matrx[i,k]=np.max(b)
        index_max_corr[i,k]=np.where(b==np.max(b))[0]

stim_matrix=np.reshape(corr_matrx,(len(corr_matrx)**2,1))
stim_matrix=stim_matrix[stim_matrix<1]

if len(names)>=2:
    plt.figure()
    plt.title('stim')
    sns.heatmap(corr_matrx, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\stim_corr.svg')

    plt.figure()
    plt.title('stim Lag')
    sns.heatmap(index_max_corr, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\stim_corr_lag.svg')

#poststim1

for i in range(len(active_indices)):
    d[header[active_indices[i]]]=a[header[active_indices[i]]]

start=int(intervall[3][2])
post=d[start:start+100]
names=list(d.columns)
names=np.array(names)[activity_overview[3,:]==1]
corr_matrx=np.zeros((len(names),len(names)))
index_max_corr=np.zeros((len(names),len(names)))
b=np.zeros((40,1))
for i in range(len(names)):
    for k in range(len(names)):
        for n in range(40):
            trace_1=(pre[names[i]])[0:stimulus_length]
            trace_2=(pre[names[k]])[0:stimulus_length]
            b[n]=crosscorr(post[names[i]],post[names[k]], lag=n)
        corr_matrx[i,k]=np.max(b)
        index_max_corr[i,k]=np.where(b==np.max(b))[0]

post_matrix_1=np.reshape(corr_matrx,(len(corr_matrx)**2,1))
post_matrix_1=post_matrix_1[post_matrix_1<1]

if len(names)>=2:
    plt.figure()
    plt.title('poststim1')
    sns.heatmap(corr_matrx, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\poststim_1_corr.svg')

    plt.figure()
    plt.title('poststim1 Lag')
    sns.heatmap(index_max_corr, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\posttim_1_corr_lag.svg')

#poststim2

for i in range(len(active_indices)):
    d[header[active_indices[i]]]=a[header[active_indices[i]]]

start=int(intervall[3][2])

post=d[-100:]
names=list(d.columns)
names=np.array(names)[activity_overview[4,:]==1]
corr_matrx=np.zeros((len(names),len(names)))
index_max_corr=np.zeros((len(names),len(names)))
b=np.zeros((40,1))
for i in range(len(names)):
    for k in range(len(names)):
        for n in range(40):
            trace_1=(post[names[i]])[0:stimulus_length]
            trace_2=(post[names[k]])[0:stimulus_length]
            b[n]=crosscorr(post[names[i]],post[names[k]], lag=n)
            corr_matrx[i,k]=np.max(b)
            index_max_corr[i,k]=np.where(b==np.max(b))[0]

post_matrix_2=np.reshape(corr_matrx,(len(corr_matrx)**2,1))
post_matrix_2=post_matrix_2[post_matrix_2<1]

if len(names)>=2:
    plt.figure()
    plt.title('poststim2')
    sns.heatmap(corr_matrx, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\poststim_2_corr.svg')

    plt.figure()
    plt.title('poststim2 Lag')
    sns.heatmap(index_max_corr, cmap='Greys')
    plt.ylim(len(corr_matrx)+0.5,0-0.5)
    plt.xlim(0-0.5,len(corr_matrx)+0.5)
    plt.ylabel('Cell ID')
    plt.xlabel('Cell ID')

    plt.savefig(folder+'\\plots\\posttim_2_corr_lag.svg')


#control

files=str(control+'\\results\\manual\\time_series.xlsx')
meta=str(control+'\\Intervals.txt')

control_df =pd.read_excel(files,sheet_name=None)


control_a=control_df['detrended']

header=list(control_a.columns.values)

v=0
peaks_all=np.zeros((len(header),1))
for i in header:
        
        zscore=stats.zscore(control_a[i])
        peaks=np.where(zscore>=2.)
        peaks=peaks[0]
        peak_diff=np.diff(peaks)
        peak_diff=peak_diff
        number_of_peaks=len(peak_diff[peak_diff>=10])
        peaks_all[v]=number_of_peaks+1
        plt.figure()
        plt.title(files[30:]+' '+i+ 'peaks='+str(number_of_peaks+1))
        plt.axis('off')
        ax=plt.subplot()
        ax.plot(control_a[i], 'k')
        ax.plot(control_a[i][peaks],'r.')
        plt.savefig(folder+'\\plots\\control_single traces\\' +i+'.svg')
        v=v+1

control=control_a[0:200]
names=header
corr_matrx=np.zeros((len(names),len(names)))
index_max_corr=np.zeros((len(names),len(names)))
b=np.zeros((40,1))
for i in range(len(names)):
    for k in range(len(names)):
        for n in range(40):
            trace_1=(control[names[i]])[0:stimulus_length]
            trace_2=(control[names[k]])[0:stimulus_length]
            b[n]=crosscorr(control[names[i]],control[names[k]], lag=n)
        corr_matrx[i,k]=np.max(b)
        index_max_corr[i,k]=np.where(b==np.max(b))[0]

control_matrix=np.reshape(corr_matrx,(len(corr_matrx)**2,1))
control_matrix=control_matrix[control_matrix<1]

plt.figure()
plt.title('control')
sns.heatmap(corr_matrx, cmap='Greys')
plt.ylim(len(corr_matrx)+0.5,0-0.5)
plt.xlim(0-0.5,len(corr_matrx)+0.5)
plt.ylabel('Cell ID')
plt.xlabel('Cell ID')

plt.savefig(folder+'\\plots\\control_corr.svg')

plt.figure()
plt.title('Control Lag')
sns.heatmap(index_max_corr, cmap='Greys')
plt.ylim(len(corr_matrx)+0.5,0-0.5)
plt.xlim(0-0.5,len(corr_matrx)+0.5)
plt.ylabel('Cell ID')
plt.xlabel('Cell ID')

plt.savefig(folder+'\\plots\\control_corr_lag.svg')


#autocorrelation pre1
start=int(intervall[3][1])
a=d[0:start]
names=list(d.columns)
corr_matrix_pre1=np.zeros((80,len(names)))
b=np.zeros((80,1))
for i in range(len(names)):
        for n in range(80):
            k=n-40
            b[n]=crosscorr(a[names[i]],a[names[i]], lag=k)
            b=np.reshape(b,(80))
            corr_matrix_pre1[:,i]=b

for i in range(len(names)):
    corr_matrx=corr_matrx[corr_matrx<1]
    plt.figure()
    plt.plot(np.linspace(-40,40,80),corr_matrix_pre1[:,i],'k-')
    sns.despine()
    #plt.xticks([])
    plt.ylabel('Autocorrelation')
    plt.xlabel('Lag [Frame]')   
    plt.savefig(folder+'\\plots\\auto\\autocorr_pre1_'+str(names[i])+'.svg')
    
    
#autocorrelation pre2
start=int(intervall[3][1])
a=d[100:200]
names=list(d.columns)
corr_matrix_pre2=np.zeros((80,len(names)))
b=np.zeros((80,1))
for i in range(len(names)):
        for n in range(80):
            k=n-40
            b[n]=crosscorr(a[names[i]],a[names[i]], lag=k)
            b=np.reshape(b,(80))
            corr_matrix_pre2[:,i]=b

for i in range(len(names)):
    #corr_matrx=corr_matrx[corr_matrx<1]
    plt.figure()
    plt.plot(np.linspace(-40,40,80),corr_matrix_pre2[:,i],'k-')
    sns.despine()
    #plt.xticks([])
    plt.ylabel('Autocorrelation')
    plt.xlabel('Lag [Frame]')   
    plt.savefig(folder+'\\plots\\auto\\autocorr_pre2_'+str(names[i])+'.svg')

#autocorrelation stim   
start=int(intervall[2][2])
stop=int(intervall[2][3])
a=d[start:stop]
names=list(d.columns)
corr_matrix_stim=np.zeros((80,len(names)))
b=np.zeros((80,1))
for i in range(len(names)):
        for n in range(80):
            k=n-40
            b[n]=crosscorr(a[names[i]],a[names[i]], lag=k)
            b=np.reshape(b,(80))
            corr_matrix_stim[:,i]=b

for i in range(len(names)):
    #corr_matrx=corr_matrx[corr_matrx<1]
    plt.figure()
    plt.plot(np.linspace(-40,40,80),corr_matrix_stim[:,i],'k-')
    sns.despine()
    #plt.xticks([])
    plt.ylabel('Autocorrelation')
    plt.xlabel('Lag [Frame]')   
    plt.savefig(folder+'\\plots\\auto\\autocorr_stim_'+str(names[i])+'.svg')
    
#autocorrelation poststim1   
start=int(intervall[3][2])
a=d[start:]
names=list(d.columns)
corr_matrix_post1=np.zeros((80,len(names)))
b=np.zeros((80,1))
for i in range(len(names)):
        for n in range(80):
            k=n-40
            b[n]=crosscorr(a[names[i]],a[names[i]], lag=k)
            b=np.reshape(b,(80))
            corr_matrix_post1[:,i]=b

for i in range(len(names)):
    #corr_matrx=corr_matrx[corr_matrx<1]
    plt.figure()
    plt.plot(np.linspace(-40,40,80),corr_matrix_post1[:,i],'k-')
    sns.despine()
    #plt.xticks([])
    plt.ylabel('Autocorrelation')
    plt.xlabel('Lag [Frame]')   
    plt.savefig(folder+'\\plots\\auto\\autocorr_poststim1_'+str(names[i])+'.svg')
    
#autocorrelation poststim2   
start=int(intervall[3][2])
a=d[-100:]
names=list(d.columns)
corr_matrix_post2=np.zeros((80,len(names)))
b=np.zeros((80,1))
for i in range(len(names)):
        for n in range(80):
            k=n-40
            b[n]=crosscorr(a[names[i]],a[names[i]], lag=k)
            b=np.reshape(b,(80))
            corr_matrix_post2[:,i]=b

for i in range(len(names)):
    #corr_matrx=corr_matrx[corr_matrx<1]
    plt.figure()
    plt.plot(np.linspace(-40,40,80),corr_matrix_post2[:,i], 'k-')
    sns.despine()
    #plt.xticks([])
    plt.ylabel('Autocorrelation')
    plt.xlabel('Lag [Frame]')
    plt.savefig(folder+'\\plots\\auto\\autocorr_poststim2_'+str(names[i])+'.svg')

#autocorr ctrl  
start=int(intervall[3][2])
a=control_a[-100:]
names=list(a.columns)
corr_matrix_autoctrl=np.zeros((80,len(names)))
b=np.zeros((80,1))
for i in range(len(names)):
        for n in range(80):
            k=n-40
            b[n]=crosscorr(a[names[i]],a[names[i]], lag=k)
            b=np.reshape(b,(80))
            corr_matrix_autoctrl[:,i]=b

for i in range(len(names)):
    #corr_matrx=corr_matrx[corr_matrx<1]
    plt.figure()
    plt.plot(np.linspace(-40,40,80),corr_matrix_autoctrl[:,i], 'k-')
    sns.despine()
    #plt.xticks([])
    plt.ylabel('Autocorrelation')
    plt.xlabel('Lag [Frame]')
       
    plt.savefig(folder+'\\plots\\auto\\autocorr_ctrl_'+str(names[i])+'.svg')


#sns.set_palette("Pastel1")
#sns.set_palette("Pastel2")
sns.set_palette("rocket_r")
#sns.set_palette("Reds")
f, (ax_swarm, ax_hist) = plt.subplots(1,2, sharey=True, gridspec_kw={"width_ratios": (.85, .15)})

plot_data=pre_matrix_1[pre_matrix_1<0.99],pre_matrix_2[pre_matrix_2<0.99],stim_matrix[stim_matrix<0.99],post_matrix_1[post_matrix_1<0.99],post_matrix_2[post_matrix_2<0.99],control_matrix[control_matrix<0.99]
plt.figure(98)
sns.violinplot(data=plot_data,ax=ax_swarm).set_title('Crosscorrelation')
#plt.sca(f[ax_swarm])
#plt.xticks((0,1,2,3,4,5), ['Pre1','Pre2','Stim','Post1','Post2', 'Ctrl'])
ax_swarm.set_xticklabels(['Pre1','Pre2','Stim','Post1','Post2', 'Ctrl'])
sns.distplot(plot_data[0],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[1],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[2],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[3],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[4],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[5],ax=ax_hist,hist=False,vertical=True)

ax_swarm.tick_params(right=False,labelright=True)
ax_swarm.tick_params(left=False,labelleft=False)
sns.despine(ax=ax_swarm, top=True,left=True,right=True)
sns.despine(ax=ax_hist, bottom=True, top=True, right=True)


print(stats.kruskal(pre_matrix_1,pre_matrix_2,stim_matrix,post_matrix_1,post_matrix_2,control_matrix))


stats_matrix_corr=np.zeros((len(plot_data),len(plot_data)))
for u in range(len(plot_data)):
    for s in range(len(plot_data)):
        stats_matrix_corr[u][s]=stats.ttest_ind(plot_data[u],plot_data[s])[1]
    



plt.savefig(folder+'\\plots\\overview.svg')

pre_matrix_1=np.reshape(corr_matrix_pre1,(np.shape(corr_matrix_pre1)[0]*np.shape(corr_matrix_pre1)[1],1))
pre_matrix_2=np.reshape(corr_matrix_pre2,(np.shape(corr_matrix_pre2)[0]*np.shape(corr_matrix_pre2)[1],1))
stim_matrix=np.reshape(corr_matrix_stim,(np.shape(corr_matrix_stim)[0]*np.shape(corr_matrix_stim)[1],1))
post_matrix_1=np.reshape(corr_matrix_post1,(np.shape(corr_matrix_post1)[0]*np.shape(corr_matrix_post1)[1],1))
post_matrix_2=np.reshape(corr_matrix_post2,(np.shape(corr_matrix_post2)[0]*np.shape(corr_matrix_post2)[1],1))
control_matrix=np.reshape(corr_matrix_autoctrl,(np.shape( corr_matrix_autoctrl)[0]*np.shape(corr_matrix_autoctrl)[1],1))

f, (ax_swarm, ax_hist) = plt.subplots(1,2, sharey=True, gridspec_kw={"width_ratios": (.85, .15)})

plot_data=pre_matrix_1[pre_matrix_1<0.99],pre_matrix_2[pre_matrix_2<0.99],stim_matrix[stim_matrix<0.99],post_matrix_1[post_matrix_1<0.99],post_matrix_2[post_matrix_2<0.99],control_matrix[control_matrix<0.99]

plt.figure(99)
sns.violinplot(data=plot_data,ax=ax_swarm).set_title('Autocorrelation')
#plt.sca(f[ax_swarm])
#plt.xticks((0,1,2,3,4,5), ['Pre1','Pre2','Stim','Post1','Post2', 'Ctrl'])

ax_swarm.set_xticklabels(['Pre1','Pre2','Stim','Post1','Post2', 'Ctrl'])
sns.distplot(plot_data[0],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[1],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[2],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[3],ax=ax_hist,hist=False,vertical=True)
sns.distplot(plot_data[4],ax=ax_hist,hist=False,vertical=True)

#if intervall[3][2]+200<intervall[3][3]:
#    sns.distplot(plot_data[5],ax=ax_hist,hist=False,vertical=True)
ax_swarm.tick_params(right=False,labelright=True)
ax_swarm.tick_params(left=False,labelleft=False)
sns.despine(ax=ax_swarm, top=True,left=True,right=True)
sns.despine(ax=ax_hist, bottom=True, top=True, right=True)


print(stats.kruskal(pre_matrix_1,pre_matrix_2,stim_matrix,post_matrix_1,post_matrix_2,control_matrix))

    
stats_matrix_auto=np.zeros((len(plot_data),len(plot_data)))
for u in range(len(plot_data)):
    for s in range(len(plot_data)):
        stats_matrix_auto[u][s]=stats.ttest_ind(plot_data[u],plot_data[s])[1]

plt.savefig(folder+'\\plots\\autooverview.svg')

stats_corr_excel=pd.DataFrame(stats_matrix_corr, columns=['PreStim1', 'Prestim2','Stim','Poststim1','Poststim2', 'Control'],index=['PreStim1', 'Prestim2','Stim','Poststim1','Poststim2', 'Control'])
stats_corr_excel.to_excel(folder+'\\stats_corr.xlsx')

stats_auto_excel=pd.DataFrame(stats_matrix_auto,  columns=['PreStim1', 'Prestim2','Stim','Poststim1','Poststim2', 'Control'],index=['PreStim1', 'Prestim2','Stim','Poststim1','Poststim2', 'Control'])
stats_auto_excel.to_excel(folder+'\\stats_auto.xlsx')

activity_overview_excel=pd.DataFrame(activity_overview, columns=list(d.columns), index=['PreStim1', 'Prestim2','Stim','Poststim1','Poststim2'])
activity_overview_excel.to_excel(folder+'\\activity.xlsx')
