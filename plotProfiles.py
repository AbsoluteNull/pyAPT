from __future__ import division
from pylab import plot, ylim, xlim, show, xlabel, ylabel, grid
from numpy import linspace, loadtxt, ones, convolve
import numpy as np
import csv
import matplotlib.pyplot as plt
import re
import os,sys
from optparse import OptionParser

import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['lines.linewidth'] = 2.0
import matplotlib.font_manager as font_manager
ticks_font = font_manager.FontProperties(size=20)


print('################################################################')
print('                     START                                      ')
print('################################################################')

	



parser = OptionParser(usage   = "usage: %prog [options] filename", \
                      version = "%prog 1.0")
parser.add_option('-i', '--il', type="float",
                  help="interface location",
                  dest="IL", default=float(75.00))


(options, filenames) = parser.parse_args()

print('######## Reading data from file ##############')

with open(filenames[0], 'rb') as f:
    reader = csv.reader(f)
    allList = list(reader)

print('######## Data read ##############')

colors     = ['m','k','#FFD700','#FF6D00','r','b']
atomLabels = ['In','Cd','S','Ga','Se','Cu']
#atomLables = ['In','O','Cd','S','Ga','Se','Cu','Na','H']


IL = options.IL

#x  = np.array([float(d) for d in allList[0] ])
x = list((np.array(range(np.shape(allList)[1]))- IL)*0.4)

fig,(ax1)=plt.subplots(1,1)

# the format of data in CSV file
# 0 -- > dist
# 1 ---> indium 
#        oxygen  
#        cadmium 
#        sulphur 
#        galium  
#        selenium
#        copper  
#        sodium  
#        hydrogen

y  = np.array([float(d)*100.0 for d in allList[1] ])
ye = np.array([float(d)*100.0 for d in allList[1+9]])
ax1.errorbar(x, y,yerr=ye,color=colors[0],label=atomLabels[0])
y  = np.array([float(d)*100.0 for d in allList[3] ])
ye = np.array([float(d)*100.0 for d in allList[3+9]])
ax1.errorbar(x, y,yerr=ye,color=colors[1],label=atomLabels[1])
y  = np.array([float(d)*100.0 for d in allList[4] ])
ye = np.array([float(d)*100.0 for d in allList[4+9]])
ax1.errorbar(x, y,yerr=ye,color=colors[2],label=atomLabels[2])
y  = np.array([float(d)*100.0 for d in allList[5] ])
ye = np.array([float(d)*100.0 for d in allList[5+9]])
ax1.errorbar(x, y,yerr=ye,color=colors[3],label=atomLabels[3])
y  = np.array([float(d)*100.0 for d in allList[6] ])
ye = np.array([float(d)*100.0 for d in allList[6+9]])
ax1.errorbar(x, y,yerr=ye,color=colors[4],label=atomLabels[4])
y  = np.array([float(d)*100.0 for d in allList[7] ])
ye = np.array([float(d)*100.0 for d in allList[7+9]])
ax1.errorbar(x, y,yerr=ye,color=colors[5],label=atomLabels[5])


#for i in range(len(atomLabels)):
#  yy  = np.array([float(d) for d in allList[i+1] ])
#  yye = np.array([float(d) for d in allList[i+1+9]])
#  ax1.errorbar(xx, yy,yerr=yye,color=colors[i],label=atomLabels[i])


ax1.set_xlabel('Distance from the interface [nm]', fontsize=20)
ax1.set_ylabel('Atomic Percentage', fontsize=20)
ax1.axhline(linewidth=2, color="k")
ax1.axvline(linewidth=1.5, linestyle='dashed',color="k")
for label in ax1.get_xticklabels():
   label.set_fontproperties(ticks_font)
for label in ax1.get_yticklabels():
   label.set_fontproperties(ticks_font)

ax1.xaxis.set_tick_params(size=8,width=2.0)
ax1.yaxis.set_tick_params(size=8,width=2.0)

#ax1.set_ylim([0,57])
#ax1.set_xlim([-30,30])
#
#
# get handles
handles, labels = ax1.get_legend_handles_labels()
## remove the errorbars
handles = [h[0] for h in handles]
# use them in the legend
ax1.legend(handles, labels, loc='upper left',numpoints=1)

###-----Plot concentration profiles -------------------
#plt.legend()
plt.show()
fig.savefig(filenames[0].split('.')[0]+'-CPall.svg',dpi=fig.dpi)
fig.savefig(filenames[0].split('.')[0]+'-CPall.pdf',dpi=fig.dpi)
fig.savefig(filenames[0].split('.')[0]+'-CPall.png',dpi=fig.dpi)
#plt.show()


print('################################################################')
print('                     END                                        ')
print('################################################################')


