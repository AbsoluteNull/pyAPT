from __future__ import division
from pylab import plot, ylim, xlim, show, xlabel, ylabel, grid
from numpy import linspace, loadtxt, ones, convolve
import numpy as np
import csv
import matplotlib.pyplot as plt
import re
import os,sys
from optparse import OptionParser
#import statsmodels.api as sm

# function for moving average calculation
def maAK(inputData, windowSize):
    padLen = int(0.01*(len(inputData)))
    p_inputData = np.lib.pad(inputData, (padLen,padLen), 'edge')
    window= np.ones(int(windowSize))/float(windowSize)
    return np.convolve(p_inputData, window, 'same')[padLen:-padLen]

def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)
#   ax.xlim([-10,10])
#   plt.show()




print('################################################################')
print('                     START                                      ')
print('################################################################')
scriptName = os.path.splitext(os.path.basename(__file__))[0]

parser = OptionParser(usage   = "usage: %prog [options] filename", \
                      version = "%prog 1.0")
parser.add_option("-m", "--movavg", type="float",
                  help="the interval width (bin size)",
                  dest="binSize", default=float(4))
parser.add_option('-i', '--il', type="float",
                  help="interface location",
                  dest="IL", default=float(20.00))
parser.add_option('-l', '--lr', type="int",
                  help="large range stuff",
                  dest="LR", default=1)
parser.add_option('-n', '--nma', type="int",
                  help="set to 1 for doing calculations without moving average",
                  dest="NMA", default=0)
(options, filenames) = parser.parse_args()
scriptName = os.path.splitext(os.path.basename(__file__))[0]
#
# syntax 'rule' for loading a file in a particular format 
print '*** Reading the data from the file',filenames[0]

(options, filenames) = parser.parse_args()
#
# syntax 'rule' for loading a file in a particular format 
with open(filenames[0], 'rb') as f:
    reader = csv.reader(f)
    allList = list(reader)

allData = []
#transpose all list
T_allList = map(list,zip(*allList))
totColumns =  len(T_allList)
 

ionNames = []
dataWbg = []

for i, row in enumerate(T_allList):
  nameInfoAll =  row[0].split()
  nameInfo    =  row[0].split()[0]
  # ignore the first few columns distance, ion cont , atom count 
  if('Sigma' not in nameInfoAll ):
  # convert the data (not the first entry that is a name)into float 
    floatData = [float(j) for j in row[1:]]
  #insert the name at the start of the list
    floatData.insert(0,nameInfo)
    dataWbg.append(floatData)

# smoothing the ion counts with moving average and inserting this column
d = dataWbg[0][1:]
icl = dataWbg[1][1:]
if (options.NMA==1):
  icnl = icl
else:
  icnl = maAK(icl,options.binSize)
  icnl = list(icnl)
icnl.insert(0,'IonS')
dataWbg.insert(2,icnl)

#
dataWObg = dataWbg[0:3]
for i, row in enumerate(dataWbg):
  nameInfoAll =  row[0].split()
  nameInfo    =  row[0].split()[0]
#  # ignore the first few columns distance, ion cont , atom count 
  if((nameInfo not in ['Distance','Ion','Atom']) and \
      'bg' not in nameInfo ):
    orig = row[1:]
    BG   = dataWbg[i+1][1:]
#    #BGsub = orig-BG 
    BGsub =  [max(0.01*c*(a-b),float(0)) for a,b,c in zip(orig,BG,icnl[1:])]
    BGsub.insert(0,nameInfo)
    dataWObg.append(BGsub)
#
overlapStuff = []
nonoverlapStuff = []

# BG substracted data and converted in counts is ready
for i, row in enumerate(dataWObg):
  nameInfoAll =  row[0].split()
  nameInfo =  row[0].split()[0]
  #ignore first few columns (alias rows)
  if((nameInfo not in ['Distance','Ion','Atom','IonS'])):
    # if encountered an overlap
    if('and' in nameInfo):
      OL  = nameInfo.split("and")
      OL2 =  OL[1].split("At")
      OL3 = OL2[1].split("bg")
      overlapSpecies = [OL[0],OL2[0]]
      overlapAt = int(OL3[0])
      overlapStuff.append([overlapSpecies,[overlapAt],row[1:]])
    else:
      nonoverlapStuff.append(row)
#
## overlap dictionaries

if (options.LR==1):
  oDict = {56:[['Cdp2','Inp2'],[42.48,4.47],[296.77,100.0],[7.27,0],[296.77,100.0]],\
         65:[['Cu','S2'],[44.78,1.6],[100.,108.96],[7.27,0],[296.77,100.0]],\
         113:[['Cd','In'],[42.48,4.47],[296.77,100.0],[7.27,0.0],[296.77,100.0]],\
         137:[['CuSe'],[1.44],[61.26],[0.0,0],[296.77,100.0]],\
         138:[['CdS'],[4.06],[211.85],[0.0,0],[296.77,100.0]],\
         140:[['CdS','CuSe'],[3.11,12.56],[211.85,61.26],[0.0,0],[296.77,100.0]],\
         142:[['CdS','CuSe'],[41.39,5.62],[211.85,61.26],[0.0,0],[296.77,100.0]],\
         143:[['CdS','CuSe'],[42.79,100.0],[211.85,61.26],[0.0,0],[296.77,100.0]],\
         145:[['CdS','CuSe'],[43.35,52.17],[211.85,61.26],[0.0,0],[296.77,100.0]],\
         147:[['CdS','CuSe'],[2.59,6.81],[211.85,61.26],[0.0,0],[296.77,100.0]]\
		}

  decoAll = []
  #
  forPlot = []
  for key, value in oDict.iteritems():
    forPlot1 = []
    for abc in overlapStuff:
      if(abc[1][0]==key):
        OVL = abc[2]
	forPlot1.append(OVL)
    for vi,v in enumerate(value[0]):
      alpha = value[1][vi]
      alpha1 = value[3][vi]
      beta = value[2][vi]
      beta1 = value[4][vi]
      for row in nonoverlapStuff:
        if(row[0]==v):
          buff = row[1:]
      deco =  [(a*alpha/beta + a*2*alpha1/beta1) for a in buff]
      forPlot1.append(deco)
      deco.insert(0,v.split('p2')[0])
      decoAll.append(deco)
    forPlot.append(forPlot1)
  #

else:
  oDict = {56:[['Cdp2','Inp2'],[42.48,4.47],[304.04,100.0]],\
         65:[['Cu','S2'],[44.78,1.6],[100.,108.96]],\
         113:[['Cd','In'],[42.48,4.47],[304.04,100.0]],\
         137:[['CuSe'],[1.44],[61.26]],\
         138:[['CdS'],[4.06],[211.85]],\
         140:[['CdS','CuSe'],[3.11,12.56],[211.85,61.26]],\
         142:[['CdS','CuSe'],[41.39,5.62],[211.85,61.26]],\
         143:[['CdS','CuSe'],[42.79,100.0],[211.85,61.26]],\
         145:[['CdS','CuSe'],[43.35,52.17],[211.85,61.26]],\
         147:[['CdS','CuSe'],[2.59,6.81],[211.85,61.26]]\
		}
  decoAll = []
  #
  for key, value in oDict.iteritems():
    for vi,v in enumerate(value[0]):
      alpha = value[1][vi]
      beta = value[2][vi]
      for row in nonoverlapStuff:
        if(row[0]==v):
          buff = row[1:]
      deco =  [(a*alpha/beta) for a in buff]
      deco.insert(0,v.split('p2')[0])
      decoAll.append(deco)
#
LOC =  len(decoAll[0])-1

indium = np.zeros(LOC)
oxygen = np.zeros(LOC)
cadmium = np.zeros(LOC)
sulphur = np.zeros(LOC)
galium = np.zeros(LOC)
selenium = np.zeros(LOC)
copper = np.zeros(LOC)
sodium = np.zeros(LOC)
hydrogen = np.zeros(LOC)
#
error_indium = np.zeros(LOC)
error_oxygen = np.zeros(LOC)
error_cadmium = np.zeros(LOC)
error_sulphur = np.zeros(LOC)
error_galium = np.zeros(LOC)
error_selenium = np.zeros(LOC)
error_copper = np.zeros(LOC)
error_sodium = np.zeros(LOC)
error_hydrogen = np.zeros(LOC)
#

for stuff in nonoverlapStuff:
  nameStuff = stuff[0]
  # decompose the name to get the constituents
  pos = [i for i,chars in enumerate(nameStuff) if chars.isupper()]
  theData = np.array(stuff[1:])
  parts=[]
  atoms = []
  stoic = []
  for j in xrange(len(pos)):
    try:
      parts.append(nameStuff[pos[j]:pos[j+1]])
    except IndexError:
      parts.append(nameStuff[pos[j]:])
  for pa in parts:
	  if(pa[-1].isdigit()):
	    if(pa[-2] is not 'p'):
              stoic.append(int(pa[-1]))
	      atoms.append(pa[:-1])
	    else:
             if(pa[-3].isdigit()):
	       stoic.append(int(pa[-3]))
	       atoms.append(pa[:-3])
             else:
	       stoic.append(1)
	       atoms.append(pa[:-2])
          else:
            stoic.append(1)
	    atoms.append(pa)
  for jj, ats in enumerate(atoms):
    if(ats=='In'):
      indium += theData*stoic[jj]
    if(ats=='Cu'):
      copper += theData*stoic[jj]
    if(ats=='Se'):
      selenium += theData*stoic[jj]
    if(ats=='Cd'):
      cadmium += theData*stoic[jj]
    if(ats=='Na'):
      sodium += theData*stoic[jj]
    if(ats=='O'):
      oxygen += theData*stoic[jj]
    if(ats=='S'):
      sulphur += theData*stoic[jj]
    if(ats=='Ga'):
      galium += theData*stoic[jj]
    if(ats=='H'):
      hydrogen += theData*stoic[jj]

for stuff in decoAll:
  nameStuff = stuff[0]
  # decompose the name to get the constituents
  pos = [i for i,chars in enumerate(nameStuff) if chars.isupper()]
  theData = np.array(stuff[1:])
  parts=[]
  atoms = []
  stoic = []
  for j in xrange(len(pos)):
    try:
      parts.append(nameStuff[pos[j]:pos[j+1]])
    except IndexError:
      parts.append(nameStuff[pos[j]:])
  for pa in parts:
	  if(pa[-1].isdigit()):
	    if(pa[-2] is not 'p'):
              stoic.append(int(pa[-1]))
	      atoms.append(pa[:-1])
	    else:
             if(pa[-3].isdigit()):
	       stoic.append(int(pa[-3]))
	       atoms.append(pa[:-3])
             else:
	       stoic.append(1)
	       atoms.append(pa[:-2])
          else:
            stoic.append(1)
	    atoms.append(pa)
  for jj, ats in enumerate(atoms):
    if(ats=='In'):
      indium += theData*stoic[jj]
    if(ats=='Cu'):
      copper += theData*stoic[jj]
    if(ats=='Se'):
      selenium += theData*stoic[jj]
    if(ats=='Cd'):
      cadmium += theData*stoic[jj]
    if(ats=='Na'):
      sodium += theData*stoic[jj]
    if(ats=='O'):
      oxygen += theData*stoic[jj]
    if(ats=='S'):
      sulphur += theData*stoic[jj]
    if(ats=='Ga'):
      galium += theData*stoic[jj]
    if(ats=='H'):
      hydrogen += theData*stoic[jj]

AAA = np.vstack((indium, \
oxygen   ,\
cadmium  ,\
sulphur  ,\
galium   ,\
selenium ,\
copper   ,\
sodium   ,\
hydrogen \
	))

TC = np.sum(AAA,axis=0)
#TC = np.shape(np.sum(AAA,axis=0))
inP = AAA/TC
intLoc = options.IL

print np.shape(AAA)

testError = np.sqrt((AAA[:,:]/TC[:])*(1 -AAA[:,:]/TC[:])/TC[:])
#testError = np.sqrt((AAA[0,:]/TC[:])*(1 -AAA[0,:]/TC[:])/TC[:])



import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 2.4
mpl.rcParams['lines.linewidth'] = 2.0
import matplotlib.font_manager as font_manager
ticks_font = font_manager.FontProperties(size=32)
from matplotlib.ticker import AutoMinorLocator


#fig1 = plt.figure()
###--- Overlap Check -------------------------
ii = [65,147,145,113,56,137,138,140,142,143]
for i,xyz in enumerate(forPlot):
  fig1 = plt.figure()
  ax2=fig1.add_subplot(111)
  plt.plot(list((np.array(d)-intLoc)*-1.0)[::-1],xyz[0][::-1],'-or',markersize=8.0,label='overlap peak')
  plt.plot(list((np.array(d)-intLoc)*-1.0)[::-1],xyz[1][1:][::-1],'-ob',markersize=8.0,label='blue')
  plt.title('overlap at '+str(ii[i]))
  if(len(xyz)>2):
    plt.plot(list((np.array(d)-intLoc)*-1.0)[::-1],xyz[2][1:][::-1],'-og',markersize=8.0,label='green')
  plt.legend(loc='upper left')
  plt.xlim([-10,10])

  ax2.xaxis.set_tick_params(size=8,width=2.0)
  ax2.yaxis.set_tick_params(size=8,width=2.0)

#  ax2.axhline(linewidth=2, color="k")
  ax2.axvline(linewidth=1.5, linestyle='dashed',color="k")
  for label in ax2.get_xticklabels():
     label.set_fontproperties(ticks_font)    
  for label in ax2.get_yticklabels():
     label.set_fontproperties(ticks_font)
 

  plt.savefig(filenames[0].split('.')[0]+str(ii[i])+'-overlapCheck.svg',dpi=fig1.dpi)
  plt.savefig(filenames[0].split('.')[0]+str(ii[i])+'-overlapCheck.png',dpi=fig1.dpi)
#  plt.show()
#
###--- Overlap Check -------------------------
def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)
    ax.xlim([-10,10])
#   plt.show()

#---------------------------------------------
subs=['In','Cd','S','Ga','Se','Cu']




###--- Plot concentration profiles -------------------------
fig,(ax1)=plt.subplots(1,1)
#
##fig1 = plt.figure()
### No error bars
##plt.plot(list(np.array(d)-intLoc),inP[0,:]*100,'-m',label='In')
##plt.plot(list(np.array(d)-intLoc),inP[2,:]*100,'-k',label='Cd')
##plt.plot(list(np.array(d)-intLoc),inP[3,:]*100,'-y',label='S')
##plt.plot(list(np.array(d)-intLoc),inP[4,:]*100,'-g',label='Ga')
##plt.plot(list(np.array(d)-intLoc),inP[5,:]*100,'-r',label='Se')
##plt.plot(list(np.array(d)-intLoc),inP[6,:]*100,'-b',label='Cu')
#
#
## with error bars
ax1.errorbar(list(-1.0*(np.array(d)-intLoc))[::-1],list(inP[0,:]*100)[::-1],yerr=list(testError[0,:]*100)[::-1],color='m',label='In')
ax1.errorbar(list(-1.0*(np.array(d)-intLoc))[::-1],list(inP[2,:]*100)[::-1],yerr=list(testError[2,:]*100)[::-1],color='k',label='Cd')
ax1.errorbar(list(-1.0*(np.array(d)-intLoc))[::-1],list(inP[3,:]*100)[::-1],yerr=list(testError[3,:]*100)[::-1],color='#FFD700',label='S')
ax1.errorbar(list(-1.0*(np.array(d)-intLoc))[::-1],list(inP[4,:]*100)[::-1],yerr=list(testError[4,:]*100)[::-1],color='#FF6D00',label='Ga')
ax1.errorbar(list(-1.0*(np.array(d)-intLoc))[::-1],list(inP[5,:]*100)[::-1],yerr=list(testError[5,:]*100)[::-1],color='r',label='Se')
ax1.errorbar(list(-1.0*(np.array(d)-intLoc))[::-1],list(inP[6,:]*100)[::-1],yerr=list(testError[6,:]*100)[::-1],color='b',label='Cu')
#
## with error fills
#errorfill(list(np.array(d)-intLoc),inP[0,:]*100,yerr=testError[0,:]*100,color='m')
#errorfill(list(np.array(d)-intLoc),inP[2,:]*100,yerr=testError[2,:]*100,color='k')
#errorfill(list(np.array(d)-intLoc),inP[3,:]*100,yerr=testError[3,:]*100,color='y')
#errorfill(list(np.array(d)-intLoc),inP[4,:]*100,yerr=testError[4,:]*100,color='g')
#errorfill(list(np.array(d)-intLoc),inP[5,:]*100,yerr=testError[5,:]*100,color='r')
#errorfill(list(np.array(d)-intLoc),inP[6,:]*100,yerr=testError[6,:]*100,color='b')

# No error bars
#plt.plot(list(np.array(d)-intLoc),inP[0,:]*100,'-m',label='In')
#plt.plot(list(np.array(d)-intLoc),inP[2,:]*100,'-k',label='Cd')
#plt.plot(list(np.array(d)-intLoc),inP[3,:]*100,'-y',label='S')
#plt.plot(list(np.array(d)-intLoc),inP[4,:]*100,'-g',label='Ga')
#plt.plot(list(np.array(d)-intLoc),inP[5,:]*100,'-r',label='Se')
#plt.plot(list(np.array(d)-intLoc),inP[6,:]*100,'-b',label='Cu')



ax1.set_xlabel('Distance from the interface')
ax1.set_ylabel('Atomic Percentage')

ax1.axhline(linewidth=2, color="k")
ax1.axvline(linewidth=1.5, linestyle='dashed',color="k")
for label in ax1.get_xticklabels():
   label.set_fontproperties(ticks_font)
for label in ax1.get_yticklabels():
   label.set_fontproperties(ticks_font)

ax1.xaxis.set_tick_params(size=9,width=2.5)
ax1.yaxis.set_tick_params(size=9,width=2.5)

ax1.set_ylim([0,58])
ax1.set_xlim([-10,10])
minor_locator = AutoMinorLocator(5)
ax1.xaxis.set_minor_locator(minor_locator)
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.tick_params(which='minor', length=5, width=2, color='k')




# get handles
handles, labels = ax1.get_legend_handles_labels()
## remove the errorbars
handles = [h[0] for h in handles]
# use them in the legend
ax1.legend(handles, labels, loc='upper left',numpoints=1)


###--- Plot concentration profiles -------------------------
##plt.legend()
fig.savefig(filenames[0].split('.')[0]+'-CPall.svg',dpi=fig.dpi)
fig.savefig(filenames[0].split('.')[0]+'-CPall.png',dpi=fig.dpi)
#fig.savefig(filenames[0].split('.')[0]+'-CPall.eps',dpi=fig.dpi)
#fig.savefig(filenames[0].split('.')[0]+'-CPall.tiff',dpi=fig.dpi)
##--- Plot concentration profiles -------------------------
plt.show()

print('################################################################')
print('                     END                                        ')
print('################################################################')

