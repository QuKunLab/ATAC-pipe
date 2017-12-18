
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import seaborn



def CorrFig(data,outname):
    data.corr().to_csv(outname+'.corr',index=True,header=True,sep='\t')
    seaborn.set_context('notebook', font_scale=1.2)
    fig1 = seaborn.clustermap(data.corr(), method='average', metric='euclidean', figsize=(12,12), cmap='RdBu_r')
    plt.setp(fig1.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(fig1.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.savefig(outname+'.corr.pdf')
    return




def scatterplot(data,label, outfig, figsize=(25,20),cmap='RdYlBu_r',label_ratio=[],n1=1,n2=2):

    if n1<=data.shape[1] and n2<= data.shape[1]:
        n1,n2=n1,n2
    else:
        n1=data.shape[1]-1
        n2=data.shape[1]
        print 'n1, n2 should smaller than %d!'%(data.shape[1])
    plt.figure(num=1, figsize=figsize)
    ax = plt.subplot(111)
    x_range=data.iloc[:,(n1-1)].max()-data.iloc[:,(n1-1)].min()
    y_range=data.iloc[:,(n2-1)].max()-data.iloc[:,(n2-1)].min()
    plt.xlim([data.iloc[:,(n1-1)].min()-x_range*0.2,data.iloc[:,(n1-1)].max()+x_range*0.3])
    plt.ylim([data.iloc[:,(n2-1)].min()-y_range*0.2,data.iloc[:,(n2-1)].max()+y_range*0.3])

    if len(label_ratio)<2:
        plt.xlabel(data.columns[n1-1],fontsize=25)
        plt.ylabel(data.columns[n2-1],fontsize=25)
    else:

        plt.xlabel(data.columns[n1-1]+'('+str(round(label_ratio[n1-1],4)*100)+'%)',fontsize=25)
        plt.ylabel(data.columns[n2-1]+'('+str(round(label_ratio[n2-1],4)*100)+'%)',fontsize=25)

    plt.xticks(rotation=45)
    colormap=getattr(plt.cm,cmap)
    colorst=[colormap(i) for i in np.linspace(0,1,len(set(label)))]
    for i in range(len(set(label))):
        x=data.loc[label==list(set(label))[i]].iloc[:,n1-1]
        y=data.loc[label==list(set(label))[i]].iloc[:,n2-1]
        ax.scatter(x,y,color=colorst[i],label=list(set(label))[i])
        for j in range(len(list(x.index))):
#            ax.scatter(x,y,color=colorst[i],
#                        label=list(set(label))[i])
            ax.annotate(x.index[j],
                         xy=(x[j],y[j]),
                         xytext=(0, -10),
                         textcoords='offset points',
                         ha='center',
                         va='top')
    plt.legend(loc = 'upper right',fontsize=25)
    plt.savefig(outfig, dpi=300)
    plt.close()

def myLabel(inf,indf):
    return pd.Series.from_csv(inf,index_col=0,sep="\t")[indf.columns]


def myPCA(data,labelfile,outname,n=3):
    n=n if n < data.shape[1] else data.shape[1]
    pca=PCA(n_components=n)
    pca.fit(data)
    pca_ratio=pca.explained_variance_ratio_
    pca_data=pd.DataFrame(np.transpose(pca.components_),index=data.columns,columns=["PC"+str(x) for x in range(1,n+1)])
    pca_data.to_csv(outname+'.pca',index=True, header=True, sep='\t')

    label=myLabel(labelfile,data)
    scatterplot(pca_data,label,figsize=(8,8),label_ratio=pca_ratio,n1=1, n2=2,outfig=outname+'.PC1_PC2.pdf')
    scatterplot(pca_data,label,figsize=(8,8),label_ratio=pca_ratio,n1=1, n2=3,outfig=outname+'.PC1_PC3.pdf')
    scatterplot(pca_data,label,figsize=(8,8),label_ratio=pca_ratio,n1=2, n2=3,outfig=outname+'.PC2_PC3.pdf')
    return
