#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 11:15:27 2017

@author: zuozuqi
"""

import pandas as pd
import numpy as np
import operator
from scipy.stats import ttest_ind as ttest
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from numpy import mean, ptp, std

def calculateSig(comlist,outdir,f=0,p=1,q=1):
    Siglist=[]
    for x in comlist:
        comdf=pd.read_table(x,header=0,index_col=0,sep=',')
        Upcomdf=comdf[(comdf['log2FoldChange']>float(f))&(comdf['pval']<float(p))&(comdf['padj']<float(q))]
        Downcomdf=comdf[(comdf['log2FoldChange']<(float(f)*(-1)))&(comdf['pval']<float(p))&(comdf['padj']<float(q))]
        Upcomdf.to_csv(x[:-4]+'.Upcom',sep='\t',header=True,index=True)
        Downcomdf.to_csv(x[:-4]+'.Downcom',sep='\t',header=True,index=True)
        Siglist.append(list(Upcomdf.index))
        Siglist.append(list(Downcomdf.index))

    SigPeaks = [y for x in Siglist for y in x]
    return SigPeaks


def BLpeaks(data,label,outdir,icv=None,iptp=None,istd=None):
    clusters=list(data.groupby(label,axis=1))
    fptp=lambda x : ptp(x)
    fstd=lambda x : std(x)
    fcv=lambda x : std(x)/(mean(x)+1)
    BLPeaks=[]
    for cn in clusters:
        tptp=cn[1].apply(fptp, axis=1)
        if iptp == None: iptp=max(tptp)
        tstd=cn[1].apply(fstd, axis=1)
        if istd == None: istd=max(tstd)
        tcv=cn[1].apply(fcv, axis=1)
        if icv == None: icv=max(tcv)
        tBLPeaks=list(data[(tptp>iptp)|(tstd>istd)|(tcv>icv)].index)
        for peak in tBLPeaks:
            if peak not in BLPeaks:
                BLPeaks.append(peak)
    return BLPeaks

def rmBLPeaks(inPeaks, BLPeaks):
    return list(set(inPeaks).difference(set(BLPeaks)))

def PeakID2Data(PeakID,data):
    return data.loc[PeakID,:]

def PeakID2SE(PeakID,Peaklist):
    tempse=Peaklist.loc[PeakID,:]
    tempse.insert(3,'PeakID',tempse.index)
    return tempse
