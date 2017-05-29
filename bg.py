#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import re
import os,sys
from optparse import OptionParser

#-----------------------------------
def mySplit(s):
     head = s.rstrip('0123456789')
     tail = s[len(head):]
     return float(tail)

def getNewName(ndata):
  usefulNames = []
  for i in range(len(ndata)):
    splitted_ndata= ndata[i].split(":",1)
    unwantedThing = "Name"
    if unwantedThing in splitted_ndata: splitted_ndata.remove(unwantedThing)
    usefulNames.append(splitted_ndata)
  ionName = ''.join([y for x in usefulNames for y in x])
  return ionName+'bg'
#-----------------------------------

print('################################################################')
print('                     START                                      ')
print('################################################################')
scriptName = os.path.splitext(os.path.basename(__file__))[0]

parser = OptionParser(usage   = "usage: %prog [options] filename", \
                      version = "%prog 1.0")
parser.add_option("-s", "--shift", type="float",
                  help="shifts the peak by shiftFactor*half-width",
                  dest="shiftFactor", default=float(3))

(options, filenames) = parser.parse_args()


for name in filenames:
  print 'Looping over provided files'
  print 'Current File: ',name
  # file object 
  fold = open(name,'r')
  fnew = open('modified'+name,'w')
  
  rangeElem = []
  numIons = 0
  numRanges = 0
  ionList = []
  n_ionList = []
  orig_rangeData = []
  new_rangeData = []
  newData = []
  oldData = []
  
  lines = fold.readlines()
  for lineNo, line in enumerate(lines):
    words = line.split()
    if "Ions" in words[0]:
      numIons = int(mySplit(lines[lineNo+1].split()[0]))
      for ionsNo in range(numIons):
        ionList.append(str(lines[lineNo+2+ionsNo].split()[0].split("=",1)[1]))
    if "Ranges" in words[0]: 
      numRanges = mySplit(lines[lineNo+1].split()[0])
      for rangessNo in range(int(numRanges)):
        lineList =  lines[lineNo+2+rangessNo].split()
        numColumns =  len(lineList)
        if(len(lineList)>=4):
        # read
          lowRange = float(lineList[0].split("=",1)[1])
          upRange  = float(lineList[1])
          volData  = lineList[2]
          nameData = lineList[3:]
          nameData = nameData[0:-1]
          colData  = lineList[-1]
        # calc new ranges
          shift = options.shiftFactor*(upRange-lowRange)*0.5
          new_lowRange = lowRange- shift
          new_upRange  = upRange - shift
        # new name and update list
          newName = getNewName(nameData)
          n_ionList.append(newName)
        # store
          oldData.append([rangessNo+1,lowRange,upRange,volData,nameData,colData])
          newData.append([rangessNo+1,new_lowRange,new_upRange,volData,newName,colData])
        else:
          sys.exit('some inconsistency in the file around line '+str(lineNo+3+rangessNo)+' in the input file')

# find unique names
  totalIonList = np.unique(ionList + n_ionList)
# just a counter
  counter = 1
# start writing
  fnew.write('%6s\n'                    % ('[Ions]'))
  fnew.write('%s%d\n'                   % ('Number=',len(totalIonList)))
  for i,ion in enumerate(totalIonList):
    fnew.write('%s\n'                   % ('Ion'+str(i+1)+'='+ion))  
  fnew.write('%s\n'                     % ('[Ranges]'))
  fnew.write('%s%d\n'                   % ('Number=',numRanges*2))
  for i in range(int(numRanges)):
    # same old
    fnew.write('%s%.4f %.4f %s %s %s\n' % ('Range'+str(counter)+'=',oldData[i][1], oldData[i][2], oldData[i][3], ' '.join(oldData[i][4]), oldData[i][-1]))
    counter = counter+1
    # new data
    fnew.write('%s%.4f %.4f %s %s %s\n' % ('Range'+str(counter)+'=',newData[i][1], newData[i][2], newData[i][3], 'Name:'+(newData[i][4]), newData[i][-1]))
    counter = counter+1
# finish writing and close files
  fnew.close()
  fold.close()

print('################################################################')
print('                     END                                        ')
print('################################################################')
