#!/usr/bin/python

import re 
import pandas as pd
import os
import Queue
import sys
import getopt
import math

		

def findclippingchance(chrom):
	
	clipfile=chrom+'_findclip.csv'
	
	t0=pd.read_csv(clipfile,sep=',',header=None,dtype='int')

	t0=t0[(t0[3]==0) | (t0[3]>=10)]

	t0=t0[(t0[1]-t0[0]>10)]
	
	print chrom,'totalreads=',len(t0)
	
	allclips=list(t0[2])

	allstarts=[x+10 if allclips[i]-x>10 else x for i,x in enumerate(t0[0])]
	
	allends=[x-10 if allclips[i]-x<-10 else x for i,x in enumerate(t0[1])]
	
	clipsize=list(t0[3])
	
	quality=list(t0[4])
	
	nums=len(allstarts)
	
	nums2=len(allstarts)+len(allclips)
	
	allcordis=allstarts+allclips+allends
		
	allindex=sorted(range(len(allcordis)), key=lambda x: allcordis[x])
	
	
	
	records=[]	
	lastclip=0
	readsnum=0
	clipsnum=1

	lastnum=1	
	all_confidence=0
	clip_confidence=0
	lastallconfidence=1
	lastclipconfidence=1
	
	for i in allindex:
		
		if i<nums:
			
			readsnum=readsnum+1
			all_confidence=all_confidence+10**(quality[i]/10.0)
		
		elif i >=nums2:
			
			readsnum=readsnum-1

			if readsnum==0:
				all_confidence=0
			else:		
				all_confidence=all_confidence-10**(quality[i-nums2]/10.0)
		
		else:
			
			currentclip=allclips[i-nums]

			if currentclip==0:
				continue

			if currentclip==lastclip:
				
				clipsnum=clipsnum+1
				clip_confidence=clip_confidence+10**(quality[i-nums]/10.0)
			else:
				records.append([lastclip,clipsnum,lastnum,100*clipsnum/lastnum, math.log10(lastclipconfidence),math.log10(lastallconfidence)])
					
				clip_confidence=10**(quality[i-nums]/10.0)
				clipsnum=1
				
			lastclip=currentclip
			lastnum=readsnum
			lastallconfidence=all_confidence			
			lastclipconfidence=clip_confidence

	records.append([lastclip,clipsnum,lastnum,100*clipsnum/lastnum, math.log10(lastclipconfidence),math.log10(lastallconfidence)])

	records=records[1:]
	
	print 'found candidates', len(records)
	
	out=pd.DataFrame.from_records(records).to_csv('{:s}_findclip2.csv'.format(chrom),mode='w',header=None, index=False, sep=',')

		
	
	




def main():
	

	allchroms=['chrY', 'chr6', 'chr12', 'chr16', 'chr7', 'chr10', 'chr14', 'chr15', 'chr8', 'chr19', 'chr9', 'chr5', 'chr18', 'chrM', 'chr17', 'chr3', 'chr13', 'chr20', 'chr11', 'chr2', 'chr1', 'chr4', 'chr21','chrX']


	for chrom in allchroms:
		
		findclippingchance(chrom)
		
		
		

	
	
if __name__=='__main__':
	
	opts,args=getopt.getopt(sys.argv[1:],"f:")
	for op, value in opts:
		if op=='-f':
			inputfile=value
	
	main()



