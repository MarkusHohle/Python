# -*- coding: utf-8 -*-
"""
Created on Thu May  4 10:22:20 2023

@author: hohle
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 13:18:35 2020

@author: hohle
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics
#from poissfit import poissfit

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

  
def poissfit(data,CI,plotstatement):

    
    eps      = 2.220446049250313e-16
    t        = len(data)
    n        = np.sum(data)
    lamGuess = np.mean(data)
    
    lamMax   = 5*(lamGuess + 1)
    dlam     = lamMax/1e4
    lamRange = np.linspace(0,lamMax, int(1e4))
      
    #-----------------------------------------------------------------------------
    if n<10:
        n_fac  = statistics.math.factorial(n)
        exp    = np.exp(-lamRange*t) 
        pdflam = (((lamRange*t)**n)/n_fac)*exp
      
    else:
        pdflamlog = n*np.log(lamRange*t + eps)- \
        sum(np.log(np.linspace(1,n,n)))- \
        lamRange*t #log P(lam|data)
        pdflam    = np.exp(pdflamlog) #P(lam|data)
        
    C         = np.trapz(pdflam,lamRange)#first y, then x 
    pdflamN   = pdflam/C #normalizing pdf for conf int
    
###############################################################################
#evaluation:
    
    idx      = np.argmax(pdflamN)
    pdfmax   = np.amax(pdflamN)
    
#right tail
    d1     = 0
    intv1  = 0
    d2     = 1
    intv2  = 0
    
    sumintv = intv1 + intv2
    
    while (sumintv < CI*0.999): #I leave a margin here in order to account for numerical inacurracy
        if (idx+d1<len(lamRange)):
            
            d1     = d1 + 1
            indi   = np.arange(idx,idx+d1,1)
            intv1  = sum(pdflamN[indi])*dlam
#left tail
        if (idx+d2>0):
            
            d2      = d2 - 1
            indi    = np.arange(idx+d2,idx,1)
            intv2   = sum(pdflamN[indi])*dlam
            
        sumintv = intv1 + intv2

    idxU    = idx + d1
    lamUP   = lamRange[idxU]
    idxL    = idx+d2
    lamLO   = lamRange[idxL]
    diff_up = lamUP - lamGuess
    diff_do = lamGuess - lamLO
###############################################################################
#####the optional ploting part#################################################

    if (plotstatement==1):
        
        xtofill     = lamRange[np.arange(idxL,idxU,1)]
        ytofill     = pdflamN[np.arange(idxL,idxU,1)]
        ytofill[0]  = 0
        ytofill[-1] = 0
        
        plt.plot(lamRange,pdflamN)
        plt.xlabel('estimated rate $\lambda$')
        plt.ylabel('P($\lambda$|data)')
        plt.fill(xtofill, ytofill, facecolor = 'black', alpha = 0.2)
        plt.plot([lamLO, lamLO],[0,pdflamN[idxL]],'k-')
        plt.plot([lamUP, lamUP],[0,pdflamN[idxU]],'k-')
        plt.plot([lamGuess, lamGuess],[0,pdfmax],'k--')
        plt.title('$\lambda = %1.2f ^{+ %1.2f}_{-%1.2f}$' % (lamGuess, \
                  diff_up, diff_do))
        
        plt.savefig('poiss.pdf',orientation = 'landscape')

    return([lamLO, lamGuess, lamUP])

