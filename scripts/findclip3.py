#!/usr/bin/python

import pandas as pd
import os
import Queue
import sys
import getopt
import math

opts,args=getopt.getopt(sys.argv[1:],"f:")
for op, value in opts:
	if op=='-f':
		inputfile=value
		
		
		
def findindex(datalist, l, x):
	
	if x>max(datalist):
		
		return l
		
		
	if x<min(datalist):
		
		return 0
	
	
	index=1/2
	find=0
	
	while find==0 and index>0:
		
		
		if datalist[index]>x:
			
			index=index/2
			
		elif datalist[index+1]<x:
			
			index=index/2+index
		
		else:
			
			if abs(datalist[index] - x) <= abs(datalist[index+1] - x):
			 
				return index
			
			else:
				
				return index+1



def findrange(datalist, r, x):
	
	l=len(datalist)
	
	init=findindex(datalist, l, x)
	
	if abs(datalist[init]-x) > r:
		
		return []
	
	
	lower_boundary,higher_boundary=init,init
	
	while lower_boundary>0 and datalist[lower_boundary-1]>x-r:
		
		lower_boundary=lower_boundary-1
	
	while higher_boundary<l-1 and datalist[higher_boundary+1]<x+r:
		
		higher_boundary=higher_boundary+1
		
	
	return range(lower_boundary,higher_boundary+1)

def findclip(data, r, x):
	
	index=findrange(list(data[0]), r, x)
	
	clips=[data.iloc[i][3] for i in index if data.iloc[i][1]>10]
	
	return sum(clips)


def findcloseranges(cordis, allclips):
		
	cordis_closeleftboundry=[x-10 for x in cordis]
	
	cordis_closerightboundry=[x+10 for x in cordis]
	
	allcordinates=cordis_closeleftboundry+allclips+cordis_closerightboundry
	
	allcordinates_sortindex=sorted(range(len(allcordinates)), key=lambda x: allcordinates[x])
	
	l0=len(cordis)
	
	l1=len(cordis)+len(allclips)
	
	records=[[] for x in cordis]
	
	currentrange=[]
	
	for index in allcordinates_sortindex:
		
		if index<l0:
			
			currentrange.append(index)
			
		elif index<l1:
			
			for i in currentrange:
				
				records[i].append(index-l0)
				
		else:
			
			currentrange=[x for x in currentrange if x!=index-l1]
			


	return records[:len(records)/2], records[len(records)/2:]


def findvalues(data,allcordinates):
	
	
	clipsites=list(data[0])

	clipnum=list(data[2])
	
	cliprates=list(data[3])
	
	clipchance0=list(data[4])
	
	clipchance1=list(data[5])
	
	clipchances=[10**x for x in clipchance0]

	notclipchances=[round(10**(x-y),4) for x,y in zip(clipchance1,clipchance0)]
	
	start_index_found, end_index_found=findcloseranges(allcordinates,clipsites)
	

	return [str(sum([clipnum[i] for i in index_found0]))+":"+str(sum([cliprates[i] for i in index_found0])) if len(index_found0)>0 else 0 for index_found0 in start_index_found], [str(math.log10(sum([clipchances[i] for i in index_found0])))+':'+str(sum([notclipchances[i] for i in index_found0])) if len(index_found0)>0 else 0 for index_found0 in start_index_found], [str(sum([clipnum[i] for i in index_found0]))+":"+str(sum([cliprates[i] for i in index_found0])) if len(index_found0)>0 else 0 for index_found0 in end_index_found] ,[str(math.log10(sum([clipchances[i] for i in index_found0])))+':'+str(sum([notclipchances[i] for i in index_found0])) if len(index_found0)>0 else 0 for index_found0 in end_index_found]

def main(inputfile):
	
	alldata={}
	for chrom in ['chrX','chrY', 'chr6', 'chr12', 'chr16', 'chr7', 'chr10', 'chr14', 'chr15', 'chr8', 'chr19', 'chr9', 'chr5', 'chr18', 'chrM', 'chr17', 'chr3', 'chr13', 'chr20', 'chr11', 'chr2', 'chr1', 'chr4', 'chr21']:
	
		data=pd.read_csv(chrom+'_findclip2.csv',sep=',',header=None)

		if chrom=='chrX':		
			alldata['Chr22']=data

		if chrom=='chrY':
			alldata['Chr23']=data

		alldata[chrom]=data
		del data
	
	SVfile=inputfile
	
	SV_data=pd.read_csv(SVfile, sep='\t')

	data=SV_data[(SV_data['type']=='Insertion') & (SV_data['#reference']!='hs38d1' )]

	data=data.loc[(data['size']>100)]

	chr_assem=list(data['#reference'])
	
	allstart_cliprate=[[] for x in xrange(len(data))]
	
	allend_cliprate=[[] for x in xrange(len(data))]
	
	allstart_clipchance=[[] for x in xrange(len(data))]
	
	allend_clipchance=[[] for x in xrange(len(data))]
	
	for chrom in alldata.keys():
		
		chrom_index=[i for i,x in enumerate(chr_assem) if x==chrom]

		if len(chrom_index)==0:
			continue
		
		chrom_data=data.iloc[chrom_index]
		
		start=list(chrom_data['ref_start'])
		
		end=list(chrom_data['ref_stop'])
	
	
		start_cliprate,start_clipchance,end_cliprate,end_clipchance,=findvalues(alldata[chrom],start+end)
	
		for i,x in enumerate(chrom_index):
			
			allstart_cliprate[x]=start_cliprate[i]
			
			allstart_clipchance[x]=start_clipchance[i]
			
			allend_cliprate[x]=end_cliprate[i]
			
			allend_clipchance[x]=end_clipchance[i]
	
	
	
	data['start_cliprate']=allstart_cliprate
	
	data['end_cliprate']=allend_cliprate
	
	data['start_clipchance']=allstart_clipchance
	
	data['end_clipchance']=allend_clipchance
	
	data.to_csv(SVfile+'_withclipping', sep='\t', index=False, mode='w' )


if __name__ == '__main__':
	
	main(inputfile)


