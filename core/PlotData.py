'''
Created  2012
@author: GieseS


Little plotting script which is called in the analysis of different mappings to an artificial reference genome.
It produces the following plots:




1) ROC Curve
2) Overview histograms for FP / TP.


'''

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import time
import pylab as p

##### HELP FUNCTIONS ####

def trapezoidal_rule(x, y):
    """Approximates the integral through the points a,b"""
    index =  [i+1 for i in xrange(len(x)-1)]
    xdiff = np.array([x[i]-x[i-1] for i in index])
    ysum = np.array([y[i]+y[i-1] for i in index])
    return(np.dot(xdiff,ysum)/2)


#### HELP FUNCTIONS END ############


""" 
abbreviations:
tp = True positives
fp = false positives
NM = number of mismtaches
mq = mapping quality
rq = readquality
subs = substitutions in artifical reference genome (ARG)
"""


            
###############################################################################            

def CalculateRoc2(dataArray,prefix,readsize,uniquehits,mappedreads,filename):
    """
    Calculates the adjusted ROC curve as well as the AUC value derived from the adjusted points
    and writes the ROC tables to .txt files. 
    """
    starttime= time.time()
    uniquehits = float(uniquehits)
    readsize = float(readsize)
    
    
    entries = len(dataArray)
        

    resultmatrix = np.arange(entries*2)
    resultmatrix = resultmatrix.reshape(2,entries)
        
    maxrq =  max(x.rq for x in dataArray)
    maxnm =  max(x.nm[0] for x in dataArray)
    maxGaps=  max(x.gaps[0] for x in dataArray)
    maxMism=  max(x.mism[0] for x in dataArray)
    
        
    minrq =  min(x.rq for x in dataArray)
    minnm =  min(x.nm[0] for x in dataArray)
    minmq=  min(x.mq[0] for x in dataArray)
    minGaps=  min(x.gaps[0] for x in dataArray)  
    minMism=  min(x.mism[0] for x in dataArray)      
    
    
    # adjust stepsize for rq since the score behaves the other way
    quants = [1,2,3,4,5]
    tempa = maxrq-minrq
    stepsize = tempa/5
        
    rqQuants = [round(minrq+(i-1)*stepsize,3) for i in quants]
    rqQuants.reverse()
    rqQuants[-1] =0 # last entry is rounded bigger than the smallest in the dataset
        
    nmQuants = [i*maxnm/5 for i in quants]
    GapsQuants = [i*maxGaps/5 for i in quants]
    MismQuants = [i*maxMism/5 for i in quants]

    rocvector = []
        
    # i = NM,l = RQ, k = MQ
    for l in quants: # RQ
                for k in quants: # GAPS
                    for j in quants: # MISMATCH
                            temparray = [m for m in dataArray if m.gaps[0] <= GapsQuants[k-1] and m.mism[0] <= MismQuants[j-1]  and m.rq >=rqQuants[l-1]]
                            

                            tempids = [m.id for m in temparray]
                            uniquereads = {}
                            for i in xrange(0,len(tempids)):
                                uniquereads[tempids[i]] = ""

                            mappedreads = len(uniquereads)
                            
 
                            
                            templength = len(temparray)
                            
                            if templength == 0:
                                continue
                            else:
                                tempTP = sum(x.mr[0] for x in temparray)
                                tempFP =templength-tempTP
                                F = round((float(mappedreads)/ readsize) ,3)
                                sens = round((tempTP/ uniquehits) * F,3)
                                if tempFP == 0:
                                    spec = 0
                                else:
                                    spec = round((tempFP / uniquehits) * F,3)                    
                                
                                rocvector.append([rqQuants[l-1],GapsQuants[k-1],MismQuants[j-1],tempTP,tempFP,templength,sens,spec,F])
                        
                    #print ("%d\t%d\t%d\t" % (templength,tempTP,tempFP))

    #0 = NM        4 = TP        7 = sens
    #1 = RQ        5 = FP        8 = 1-spec
    #2 = GAPS        6 = P        9 = F
    #append needed for last entry in AUC calculation
    rocvector.append([0,0,0,0,0,0,0,0,0])                
    nproc = np.array(rocvector)
    
    #write the sens and specificity values from nproc according to the enumeration in line 149. 
    #specificity is in cell -2
    # sensitivity is in cell -3
    sens =  [i[-3] for i in nproc]
    spez =  [i[-2] for i in nproc]
    
    # adjust ROC curve. It is necessary that it the 1-specificity ends in 1.
    # for the last record copy the  predecessor in sens to it
    # and write 1 to specificity    
    spez[-1] = 1
    sens[-1] = sens[-2]
    

    rocarray1 = np.array([sens,spez])
    rocarray1 = rocarray1.flatten('F')
    rocarray1= rocarray1.reshape((len(spez),2))
    
    rocarray = np.array([sens,spez])
    rocarray = rocarray.flatten('F')
    rocarray = rocarray.reshape((len(spez),2))
    rocarray = np.sort(rocarray.view('float,float'), order=['f0','f1'], axis=0).view(np.float)
    
    rocarrayCorrected = rocarray
    
    #print rocarrayCorrected
    # project points where...
    for m in range(len(rocarrayCorrected)-2,-1,-1):
        if (rocarrayCorrected[m,1] >= rocarrayCorrected[m+1,1]):
            rocarrayCorrected[m,1] = rocarrayCorrected[m+1,1]

            
    #print rocarrayCorrected            
    plt.hold(True)
    plt.figure()
    plt.subplot(111)
    #plt.scatter(spez, sens, c='b', marker='o', facecolor='red')
    #plt.plot(rocarray[:,1], rocarray[:,0]
    plt.plot(rocarrayCorrected[:,1],rocarrayCorrected[:,0], marker='o', markersize=7,linestyle='--', color='r', label='projected')
    plt.plot(rocarray1[:,1], rocarray1[:,0], linestyle="None",label='real',marker='.',color='g')
    plt.xlabel('1-specificity')
    plt.ylabel('sensitivity')
    plt.title(r'ROC:'+filename)
    plt.axis([-0.1,1.1,-0.1,1.1])
    plt.grid(True)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(prefix + "_ROC.pdf",format='pdf')
    plt.clf  
    
    
    AUC = trapezoidal_rule(rocarrayCorrected[:,1], rocarrayCorrected[:,0])
    
    fobj =  open(prefix+"_roctable.txt","w")
    fobj.write("RQ\tGAPS\tMM\tPTP\tFP\tP\tSn\t1-Sp\tF\r\n")
    for i in xrange(0,len(rocvector),1):
        temp = [str(k) for k in rocvector[i]]
        tempstr = "\t".join(temp)
        fobj.write(tempstr+"\r\n")

    endtime= time.time()
    return(round(AUC,3))

def plotOverviewHist(fp,tp,label,prefix,mappernames):
    """Plots true positives and false positives into 2 different histogram subplots. """
    prefix2 = "/".join(prefix.split("/")[0:-1])+"/"
    fobj = open(prefix2+"indexMappingTools.txt","w")
    for i in range(0,len(label)):
        fobj.write("%s - %s\r\n" %(i+1,mappernames[i]))
    fobj.close()
    
    x = [i for i in range(1,len(fp)*3,3)]
    xmax = max(x)+1
    ymaxTP = max(tp)+0.1
    ymaxFP = max(fp)+0.1

    ##### SUBPLOT NUMBER OF MISMATCHES ####
    y = tp
    x =x
    z = fp

 
    fig = p.figure()
    # only plot every 2nd label
    if len(label) <= 7:
        widthp = 0.7
        ticks = label
    else:
        widthp = 0.3
        ticks = [i if i%2 == 0 else "" for i in label]
    # Here we're adding 2 subplots.  The grid is set
    # up as one row, two columns.
    ax1 = fig.add_subplot(1,2,1)
    ax1.bar(x,y,width=widthp, facecolor='darkgreen')
    ax1.set_ylabel('#TP hits')
    ax1.set_xlabel('index mapping tool')
    ax1.set_title("Global comparison #TP hits")
    p.xticks(x,ticks)
    p.grid(True)
    p.axis([0,xmax,0,ymaxTP+ymaxTP*10/100])
    
    # on the second axis, make the width smaller (default is 0.8)
    ax2 = fig.add_subplot(1,2,2)
    ax2.bar(x,z,width=widthp, facecolor='darkred')
    ax2.set_ylabel('#FP hits')
    ax2.set_xlabel('index mapping tool')
    ax2.set_title("Global comparison #FP hits")
    p.axis([0,xmax,0,ymaxFP+ymaxFP*10/100])
    p.xticks(x,ticks)
    p.grid(True)
    plt.tight_layout()
    p.savefig(prefix2 + "Overall_histabs.pdf",format='pdf')
    p.clf()
    
    
    tpsum =sum(tp)
    fpsum =sum(fp)
    y = [i/float(tpsum) for i in tp]
    x =x
    z = [i/float(fpsum) for i in fp]
    ymax = max(max(z),max(y))+0.2
    
    fig = p.figure()
    # only plot every 2nd labelare provided
    if len(label) <= 7:
        ticks = label
    else:
        ticks = [i if i%2 == 0 else "" for i in label]
    # Here we're adding 2 subplots.  The grid is set
    # up as one row, two columns.
    ax1 = fig.add_subplot(1,2,1)

   
    ax1.bar(x,y,width=widthp, facecolor='darkgreen')
    ax1.set_ylabel('%TP hits')
    ax1.set_xlabel('index mapping tool')
    ax1.set_title("Global comparison %TP hits")
    p.xticks(x,ticks)
    p.grid(True)
    p.axis([0,xmax,0,1.1])
    
    # on the second axis, make the width smaller (default is 0.8)
    ax2 = fig.add_subplot(1,2,2)
    ax2.bar(x,z,width=widthp, facecolor='darkred')
    ax2.set_ylabel('%FP hits')
    ax2.set_xlabel('index mapping tool')
    ax2.set_title("Global comparison %FP hits")
    p.axis([0,xmax,0,1.1])
    p.xticks(x,ticks)
    p.grid(True)
    plt.tight_layout()
    p.savefig(prefix2 + "Overall_histper.pdf",format='pdf')
    p.clf()
    

 

