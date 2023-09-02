# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 13:18:35 2020

@author: hohle
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from poissfit import poissfit

def my_data_analysis(filename,sample_option,S):
    
    
    workbook  = pd.ExcelFile(filename)
    names     = workbook.sheet_names
    sheet1    = workbook.parse(names[0])
    col_names = sheet1.columns
    counts    = sheet1.loc[3:][col_names[9:]]
    counts    = counts.dropna()
    counts    = np.array(counts)
    
    sam       = sample_option - 1
    WT        = counts[:,0:3]
    TR        = counts[:,3 + 3*sam:6 + 3*sam]
    
    n     = len(WT)
    
    LamWT = np.zeros((n,3),float)
    LamTR = np.zeros((n,3),float)
    
    for i in range(n):
        
        vecWT = WT[i,:]
        vecTR = TR[i,:]
        
        lamWT = poissfit(vecWT,0.68,0)
        lamTR = poissfit(vecTR,0.68,0)
        
        LamWT[i,:] = lamWT
        LamTR[i,:] = lamTR
        
    
    Diff       = LamWT[:,1] - LamTR[:,1]
    
    erLamlowWT = LamWT[:,1] - LamWT[:,0] 
    erLamlowTR = LamTR[:,1] - LamTR[:,0] 
    
    erLamupWT  = LamWT[:,2] - LamWT[:,1] 
    erLamupTR  = LamTR[:,2] - LamTR[:,1]
    
    below = np.argwhere(Diff>0)
    below = below[:,0] #extracting the indices
    above = np.argwhere(Diff<=0)
    above = above[:,0]
    
        #determining sigma
    Sigma        = np.zeros((n,1),float)
    
    dbelow       = Diff[below]
    dabove       = Diff[above]
    sibelow      = np.sqrt(erLamlowWT[below]**2 + erLamupTR[below]**2)
    siabove      = np.sqrt(erLamupWT[above]**2 + erLamlowTR[above]**2)
    
    Sigma[below,0] = dbelow/sibelow
    Sigma[above,0] = -dabove/siabove
    
    significant    = np.argwhere(Sigma>=S)
    significant    = significant[:,0]
    
    Chi2_red = np.sum(Sigma**2)/(n-1)
    
    #plotting all data points
    plt.loglog()
    plt.xlabel('counts WT')                          #to avoid overlay with 
    plt.ylabel('counts TR')                          #next plot
    plt.errorbar(LamWT[:,1],LamTR[:,1],\
                        yerr = [erLamlowTR, erLamupTR],\
                        xerr = [erLamlowWT, erLamupWT],\
                        fmt='ko', alpha = 0.2)
    
    #h = plt.loglog(LamWT[significant,1],LamTR[significant,1],'ro')
    h = plt.errorbar(LamWT[significant,1],LamTR[significant,1],\
                yerr = [erLamlowTR[significant], erLamupTR[significant]],\
                xerr = [erLamlowWT[significant], erLamupWT[significant]],\
                fmt='ro')
    plt.legend(h,['significance > '+ str(S) + '$\sigma$'], frameon = False)
    plt.title('$\chi^2_{red}$ = %1.2f' % Chi2_red)
    
    plt.savefig('peptides.pdf', orientation= 'landscape')
    
    plt.show(block=False)
    
    return(LamWT,LamTR,Chi2_red)