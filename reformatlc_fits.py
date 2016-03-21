#!usr/bin/python

#This will attempt to reformat the Light curve files provided by chris to conform to the CAOM requirements
# History: vs - 09 Feb 2016
#code modified from david's reformatMOST.py script

import copy
import glob
import os
import astropy
import sys
import vos
import string
import subprocess
from astropy.io import fits
import shutil
from os import listdir
import numpy as np
from astropy.table import Table
from os.path import basename
from astropy.io import ascii

sys.path.append('source/')
import rwlcfile
import getFitsHeader

# function used to fetch FITS header information
# This is code by D Bohlender
#def getFitsHeader(url):
    #subp = subprocess.Popen(['curl', '--location-trusted', '-n', url],
            #stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #curlStdOut, curlStdErr = subp.communicate()
    #lines = curlStdOut.rstrip().split('\n')
    #header = {}
    #for line in lines:
        #keyword = line[:8].strip()
        #value = ''
        #comment = ''
        #if keyword == 'COMMENT':
            #continue
        #else:
            #try:
                #(value, comment) = line[10:].split('/')
            #except ValueError:
                #value = line[10:]
            #value = value.strip('\' ')
        #header[keyword] = value,comment

    #return header

def rwlcfile(lcfile,outfname,phdr):
  "This short fn reads the lcdata, adds keywords and write the output file. VS - 28 Feb 2016"
  
  #read the light curve data
  vC.copy(os.path.join(voMOST,'most2caomlc',lcfile),os.path.join(os.getcwd(),os.path.basename(lcfile)))
  lcdata = ascii.read(os.path.join(os.getcwd(),os.path.basename(lcfile)))

  #add keywords UTSTART and UTEND

  phdr[0].header['UTSTART'] = (lcdata[0][0],'Start date of observations')
  phdr[0].header['UTEND'] = (lcdata[len(lcdata)-1][0],'End date of observations')
 
  bhdr = fits.BinTableHDU(data=lcdata.as_array())
  phdulist = fits.HDUList([phdr[0],bhdr])
  phdulist.writeto(os.path.join(os.getcwd(),os.path.basename(outfname)),clobber='True')#write locally first
  
  #now copy it to vospace
  vC.copy(os.path.join(os.getcwd(),os.path.basename(outfname)),outfname)
  print 'File written: '+outfname
  

#cmd = 'getCert --daysValid=25 --dest=$HOME/.ssl/cadcproxy.pem'
#os.system(cmd)


vC = vos.Client()
voMOST = 'vos://cadc.nrc.ca~vospace/MOST'
voHome = os.path.join(voMOST,'most2caomlc')

#vospace data server
dws = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub'

#Collect files; Read the list of directories from the 'lightcurve_data.txt'.

#lclist = os.path.join(voHome,'resource/lightcurve_data.txt')

lclistpath = os.path.join(voHome,'resource','lightcurve_data.txt')

lclist = vC.open(lclistpath,mode='r',view='data')

#with open(lclist,'r') as f:
   #Files = sorted(f.readlines())
Files = sorted((vos.VOFile.read(lclist)).split('\n'))
#print Files[0]

#when the lc files are read and split into a list, since the last character is \n,
  #there is an empty string at the beginning. This line removes that string

Files = Files[1:]

campaignName = []
for file in Files:
  x = os.path.split(os.path.split(os.path.dirname(file))[1])[1]
  campaignName.append(x)
  

#get the corresponding target directory
targetDir, ind = np.unique(campaignName,return_index=True)
ind = np.append(ind,len(Files)-1)


#individual campaigns
for i in range(len(ind)):
  print "Processing Directory :",targetDir[i]
  destDir = os.path.join(voMOST,targetDir[i],'Processed','Lightcurves')
  
  #get information from the target log
  #log_file = vC.copy(os.path.join(voHome,'resource/target_log.csv'),os.getcwd())
  
 ### !!!!!!!!!!!!!!!#this is short term solution to keep a copy of the log file in the working directory
  log_file = os.path.join(os.getcwd(),'target_log.csv')
  
  #log_file = 'resource/target_log.csv'
  fileRunNo = np.genfromtxt(log_file,skip_header=1,delimiter=",",dtype=None,usecols=(0))

  run_no, Ptarget, yr = string.split(targetDir[i],'_')

  lines2read = (np.where(fileRunNo == run_no.lstrip("0")))[0]
  skip = len(fileRunNo)-(lines2read[0]+len(lines2read))

  names = [ 'obs_code','obs_mode','pos','primary_ID']
  Tdat = np.genfromtxt(log_file,
		     skip_header=1,
		     delimiter=",",
		     dtype=None,
		     names=names,
		     skiprows=lines2read[0]+1,
		     skip_footer=skip.astype(int),
		     usecols=(1,2,3,4))

  #Get lc files
  lcfiles = Files[ind[i]:ind[i+1]]
  
  #Process fabry image
  
  lcind = [j for j, k in enumerate(lcfiles) if 'fab' in k]
  if lcind:
    tind = (np.where(Tdat['obs_mode'] == 'fab'))[0]
  
    Tname = Tdat['primary_ID'][tind][0]
    
    hdrdir = os.path.join(targetDir[i],'Reformatted',Tname)
    hdrdir = hdrdir.translate(None,string.whitespace) #remove if there are any white spaces
    hdrdir = hdrdir.translate(None,'-')
    hdrfilename=vC.glob(os.path.join(voMOST,hdrdir,'*/*.fits'))[0]
    
    vC.copy(hdrfilename,os.path.join(os.getcwd(),os.path.basename(hdrfilename)))#easier if copied to the CWD
    phdr = fits.open(os.path.join(os.getcwd(),os.path.basename(hdrfilename)))
    
    #url = os.path.join(dws,'vospace/MOST',hdrdir,(hdrfilename.rsplit('/'))[-2],(hdrfilename.rsplit('/'))[-1])+'?fhead=true'
    #phdr = getFitsHeader(url)
    
    
    lcfile = os.path.join(voHome,lcfiles[lcind[0]])

    #depending on whether it is published or not, change the suffix
    if 'pub' in os.path.basename(lcfile):
      outfname = os.path.join(destDir,(((run_no+'_'+(Tname)+'_'+yr+'_pub_lc.fits').translate(None,string.whitespace))).translate(None,'-'))
    else:
      os.path.join(destDir,(((run_no+'_'+(Tname)+'_'+yr+'_lc.fits').translate(None,string.whitespace))).translate(None,'-'))

    rwlcfile(lcfile,outfname,phdr)
    
    #remove the fab lc file to raw folder and remove the variable from the lc list
    vC.move(lcfile,os.path.join(os.path.join(voMOST,targetDir[i],'Raw')))
    lcfiles.remove(lcfiles[lcind[0]])
    
  Tdat = np.delete(Tdat, np.where(Tdat['obs_mode'] == 'fab'))
  
  
  #Process direct images
  
  lcind = [j for j, k in enumerate(lcfiles) if 'gs' not in k]
  if lcind:
    tind = (np.where(Tdat['obs_mode'] == 'di'))[0]

    for x in tind:
      Tname = (Tdat['primary_ID'][x].translate(None,string.whitespace)).translate(None,'-')
      hdrdir = os.path.join(targetDir[i],'Reformatted',(Tname))
      hdrdir = hdrdir.translate(None,string.whitespace)
      hdrdir = hdrdir.translate(None,'-')
      hdrfilename=vC.glob(os.path.join(voMOST,hdrdir,'*/*.fits'))[0]
    
      vC.copy(hdrfilename,os.path.join(os.getcwd(),os.path.basename(hdrfilename)))#easier if copied to the CWD
      phdr = fits.open(os.path.join(os.getcwd(),os.path.basename(hdrfilename)))
      
      #depending on whether it is published or not, change the suffix
      if 'pub' in os.path.basename(lcfile):
	outfname = os.path.join(destDir,(((run_no+'_'+(Tname)+'_'+yr+'_pub_lc.fits').translate(None,string.whitespace))).translate(None,'-'))
      else:
	os.path.join(destDir,(((run_no+'_'+(Tname)+'_'+yr+'_lc.fits').translate(None,string.whitespace))).translate(None,'-'))

      rwlcfile(lcfile,outfname,phdr)
    
      #remove the fab lc file to raw folder and remove the variable from the lc list
      vC.move(lcfile,os.path.join(os.path.join(voMOST,targetDir[i],'Raw')))
      
    lcfiles = np.delete(lcfiles,tind)
    
  Tdat = np.delete(Tdat, np.where(Tdat['obs_mode'] != 'gs'))
  
  #Process Guide stars
  lcind = [j for j, k in enumerate(lcfiles) if 'gs' in k]
  
  if lcind:
    for x in range(len(Tdat)):
      Tname = (Tdat['primary_ID'][x].translate(None,string.whitespace)).translate(None,'-')
      #when there is no guide star folder take information from the log_file
      
      hdrdir = os.path.join(targetDir[i],'Reformatted','Guide_stars')
      phdr = fits.open(glob.glob(os.path.join(hdrdir,('*'+Tname+'*gs.fits')))[0])
      
      
      hdrfilename=vC.glob(os.path.join(voMOST,hdrdir,'*/*.fits'))[0]
    
      vC.copy(hdrfilename,os.path.join(os.getcwd(),os.path.basename(hdrfilename)))#easier if copied to the CWD
      phdr = fits.open(os.path.join(os.getcwd(),os.path.basename(hdrfilename)))
      
      #depending on whether it is published or not, change the suffix
      if 'pub' in os.path.basename(lcfile):
	outfname = os.path.join(destDir,(((run_no+'_'+(Tname)+'_'+yr+'_pub_lc.fits').translate(None,string.whitespace))).translate(None,'-'))
      else:
	os.path.join(destDir,(((run_no+'_'+(Tname)+'_'+yr+'_lc.fits').translate(None,string.whitespace))).translate(None,'-'))

      rwlcfile(lcfile,outfname,phdr)
    
      #remove the fab lc file to raw folder and remove the variable from the lc list
      vC.move(lcfile,os.path.join(os.path.join(voMOST,targetDir[i],'Raw')))
      
      
print 'Finished processing ',len(Files)































