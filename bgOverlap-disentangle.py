from __future__ import division
from pylab import plot, ylim, xlim, show, xlabel, ylabel, grid
from numpy import linspace, loadtxt, ones, convolve
import numpy as np
import csv
import matplotlib.pyplot as plt
import re
import os,sys
from optparse import OptionParser
import math
#import statsmodels.api as sm




################################################################
# function for moving average calculation
################################################################
def nonLocal(source,alpha,N,length):
	freq     = np.zeros(N)
	I = complex(0,1)                     # Initialize complex I
	# EVEN!!!!!!!
	for ix in range(N):
		if ix <= N/2:
			freq[ix] = float(ix)
		if ix > N/2:
			freq[ix] = float(ix-N)
	#print freq
	#Source in frequency space
	source_f = np.fft.fft(source)
	nl_f = np.zeros(N, dtype=complex)
	for i in range(len(source_f)):
		#Solve second order nonlocal smoothing
		nl_f[i]     = source_f[i] / (1.+alpha*4.*3.141**2./length**2.*abs(freq[i])**2.)

	#NL term in frequency space
	nl = np.real(np.fft.ifft(nl_f))
	return nl
################################################################



print('################################################################')
print('                     START                                      ')
print('################################################################')


################################################################
# commandline help
################################################################

parser = OptionParser(usage   = "usage: %prog [options] filename", \
                      version = "%prog 1.0")

################################################################
# setting commandline arguments
################################################################

parser.add_option("-m", "--movavg", type="float",
                  help="the interval width (bin size)",
                  dest="binSize", default=float(14))
#parser.add_option('-i', '--il', type="float",
#                  help="interface location",
#                  dest="IL", default=float(20.00))
parser.add_option('-l', '--lr', type="int",
                  help="large range stuff",
                  dest="LR", default=1)
#parser.add_option('-n', '--nma', type="int",
#                  help="set to 1 for doing calculations without moving average",
#                  dest="NMA", default=0)


(options, filenames) = parser.parse_args()
# syntax 'rule' for loading a file in a particular format 
print '*** Reading the data from the file: --> ',filenames[0]


# syntax 'rule' for loading a csv file
with open(filenames[0], 'rb') as f:
  reader = csv.reader(f)
  FD_rWise = list(reader)

FD_cWise = map(list,zip(*FD_rWise))
 
ionNames = []
dataWbg = []

################################################################
# removing the unwanted columns of original data
# row here are the columns in original data

for i, row in enumerate(FD_cWise):
  # the first entry in columnwise data is the header containing some info as strings
  nameInfoAll =  row[0].split()
  nameInfo    =  row[0].split()[0]
  # ignore the first few columns distance, ion cont , atom count 
  # and Sigma entries
  if('Sigma' not in nameInfoAll ):
    # convert the data (not the first entry that is a name) into float 
    floatData = [float(j) for j in row[1:]]
    # insert the name at the start of the list
    floatData.insert(0,nameInfo)
    dataWbg.append(floatData)

################################################################
################################################################

# getting the distance column
# getting the Ion count column and smoothing it

## smoothing the ion counts with moving average and inserting this column
##d = list(np.array(len(dataWbg[0][1:]))*0.4)
d = np.array(range(0,len(dataWbg[0][1:])))*5
icl = dataWbg[1][1:]
#icnl = nonLocal(icl,0.0,len(icl),len(icl))
icnl = icl

#icnl [:4] = icl[:4]
#icnl [-4:] = icl[-4:]

icnl.insert(0,'IonS')
dataWbg.insert(2,icnl)
#
#fig, ax1 = plt.subplots()

#ax1.plot(d,dataWbg[1:],'-or')
plt.plot(d,dataWbg[1][1:],'-or')
plt.plot(d,dataWbg[2][1:],'-b')
#plt.plot(d,icnl1,'-g')
#plt.plot(d,icnl2,'-k')
#plt.plot(d,icnl3,'-m')

plt.show()



################################################################
# now we can want to remove the background form each species
# dataWObg will hold that info
# first 4 columns will be same and then the correct species count (%) will be appended

dataWObg = dataWbg[0:3]

for i, row in enumerate(dataWbg):
  nameInfoAll =  row[0].split()
  nameInfo    =  row[0].split()[0]
  # ignore the first few columns distance, ion cont , atom count and bg columns 
  if((nameInfo not in ['Distance','Ion','Atom']) and \
      'bg'     not in nameInfo ):
    orig = row[1:]
    BG   = dataWbg[i+1][1:]
    #BGsub = orig-BG 
    BGsub =  [max(0.01*c*(a-b),float(0)) for a,b,c in zip(orig,BG,icnl[1:])]
    BGsub.insert(0,nameInfo)
    dataWObg.append(BGsub)
##

################################################################
overlapStuff = []
nonoverlapStuff = []
#
## BG substracted data and converted in counts is ready
for i, row in enumerate(dataWObg):
  nameInfoAll =  row[0].split()
  nameInfo    =  row[0].split()[0]
  #ignore first few columns (alias rows)
  if((nameInfo not in ['Distance','Ion','Atom','IonS'])):
    # if encountered an overlap
    # do some string manipulation to ........
    if('and' in nameInfo):
      OL  = nameInfo.split("and")
      OL2 =  OL[1].split("At")
      OL3 = OL2[1].split("bg")
      overlapSpecies = [OL[0],OL2[0]]
      overlapAt = int(OL3[0])
      overlapStuff.append([overlapSpecies,[overlapAt],row[1:]])
    else:
      nonoverlapStuff.append(row)
## cluster
## overlap dictionaries
#

# keys --> peakLoc
if (options.LR==1):
  oDict = { 56:[['Cdp2','Inp2'],[42.48,4.47] ,[296.77,100.0],[7.27,0]  ,[296.77,100.0]],\
            65:[['Cu','S2']    ,[44.78,1.6]  ,[100.,108.96] ,[7.27,0]  ,[296.77,100.0]],\
           113:[['Cd','In']    ,[42.48,4.47] ,[296.77,100.0],[7.27,0.0],[296.77,100.0]],\
           137:[['CuSe']       ,[1.44]       ,[61.26]       ,[0.0,0]   ,[296.77,100.0]],\
           138:[['CdS']        ,[4.06]       ,[211.85]      ,[0.0,0]   ,[296.77,100.0]],\
           140:[['CdS','CuSe'] ,[3.11,12.56] ,[211.85,61.26],[0.0,0]   ,[296.77,100.0]],\
           142:[['CdS','CuSe'] ,[41.39,5.62] ,[211.85,61.26],[0.0,0]   ,[296.77,100.0]],\
           143:[['CdS','CuSe'] ,[42.79,100.0],[211.85,61.26],[0.0,0]   ,[296.77,100.0]],\
           145:[['CdS','CuSe'] ,[43.35,52.17],[211.85,61.26],[0.0,0]   ,[296.77,100.0]],\
           147:[['CdS','CuSe'] ,[2.59,6.81]  ,[211.85,61.26],[0.0,0]   ,[296.77,100.0]]\
          }

  decoAll = []
  forPlot = []
  # loop over all the overlap locations (key s)
  for key, value in oDict.iteritems():
    #print key
    forPlot1 = []
    for abc in overlapStuff:
      if(abc[1][0]==key):
        oveLapData = abc[2] # take the orig data
	forPlot1.append(oveLapData)
    # Loop over species
    for vi,v in enumerate(value[0]):
      alpha  = value[1][vi] # constants related to natural abundance
      alpha1 = value[3][vi]
      beta   = value[2][vi]
      beta1  = value[4][vi]
      for row in nonoverlapStuff:
        if(row[0]==v):
          buff = row[1:]
      deco =  [(a*alpha/beta + a*2*alpha1/beta1) for a in buff]
      forPlot1.append(deco)
      deco.insert(0,v.split('p2')[0])
      decoAll.append(deco)
    forPlot.append(forPlot1)
  #print forPlot[0]

###
###else:
###  oDict = { 56:[['Cdp2','Inp2'],[42.48,4.47] ,[304.04,100.0]],\
###            65:[['Cu','S2']    ,[44.78,1.6]  ,[100.,108.96]] ,\
###           113:[['Cd','In']    ,[42.48,4.47] ,[304.04,100.0]],\
###           137:[['CuSe']       ,[1.44]       ,[61.26]]       ,\
###           138:[['CdS']        ,[4.06]       ,[211.85]]      ,\
###           140:[['CdS','CuSe'] ,[3.11,12.56] ,[211.85,61.26]],\
###           142:[['CdS','CuSe'] ,[41.39,5.62] ,[211.85,61.26]],\
###           143:[['CdS','CuSe'] ,[42.79,100.0],[211.85,61.26]],\
###           145:[['CdS','CuSe'] ,[43.35,52.17],[211.85,61.26]],\
###           147:[['CdS','CuSe'] ,[2.59,6.81]  ,[211.85,61.26]]\
###          }
###  decoAll = []
###  #
###  for key, value in oDict.iteritems():
###    for vi,v in enumerate(value[0]):
###      alpha = value[1][vi]
###      beta = value[2][vi]
###      for row in nonoverlapStuff:
###        if(row[0]==v):
###          buff = row[1:]
###      deco =  [(a*alpha/beta) for a in buff]
###      deco.insert(0,v.split('p2')[0])
###      decoAll.append(deco)
####


###################################################################
# Writing Overlap check data to files
###################################################################
overlapLoc = [65,147,145,113,56,137,138,140,142,143]
for i,olc_data in enumerate(forPlot):
  OLdata = np.vstack([np.array(olc_data[0][0:]),\
                      np.array(olc_data[1][1:]) \
                     ])
  if(len(olc_data)>2):
    OLdata = np.vstack([np.array(olc_data[0][0:]),\
                        np.array(olc_data[1][1:]),\
                        np.array(olc_data[2][1:])
                       ])
  olcFname = 'overlap-at-'+str(overlapLoc[i])+'.csv'
  np.savetxt(olcFname,OLdata,delimiter=",")

###################################################################
###################################################################


LOC =  len(decoAll[0])-1

indium   = np.zeros(LOC)
oxygen   = np.zeros(LOC)
cadmium  = np.zeros(LOC)
sulphur  = np.zeros(LOC)
galium   = np.zeros(LOC)
selenium = np.zeros(LOC)
copper   = np.zeros(LOC)
sodium   = np.zeros(LOC)
hydrogen = np.zeros(LOC)
#
error_indium   = np.zeros(LOC)
error_oxygen   = np.zeros(LOC)
error_cadmium  = np.zeros(LOC)
error_sulphur  = np.zeros(LOC)
error_galium   = np.zeros(LOC)
error_selenium = np.zeros(LOC)
error_copper   = np.zeros(LOC)
error_sodium   = np.zeros(LOC)
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
hydrogen ))

TC = np.sum(AAA,axis=0)
inP = AAA/TC
testError = np.sqrt((AAA[:,:]/TC[:])*(1 -AAA[:,:]/TC[:])/TC[:])
atoErr = np.vstack([d,inP,testError])
np.savetxt("atomAll.csv", atoErr, delimiter=",")

print np.shape(forPlot[0])
print('################################################################')
print('                     END                                        ')
print('################################################################')

