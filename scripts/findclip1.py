#!/usr/bin/python

import re 
import pandas as pd
import os
import Queue
import sys
import getopt
		

def readcigar(cigar,posi):

	alllabels=re.findall(r'[0-9]*[A-Z=]',cigar)
	
	if cigar.count('S')+ cigar.count('H')==1:
		
		clip=[i for i,x in enumerate(alllabels) if x[-1] in ['S','H']][0]
		
		clipsize=int(alllabels[clip][:-1])

		if clip==0:
			frontposi=posi		

		else:	
			front=alllabels[:clip]
		
			frontposi=posi+sum([int(x[:-1]) if len(x)>1 else 1 for x in front if x[-1] in ['M','I','N','X','=']])-sum([int(x[:-1]) if len(x)>1 else 1 for x in front if x[-1] in ['P','D']])		
			
	else:
		
		frontposi=0
		
		clipsize=0
		
	end=posi+sum([int(x[:-1]) if len(x)>1 else 1 for x in alllabels if x[-1] in ['M','I','N','X','=']])-sum([int(x[:-1]) if len(x)>1 else 1 for x in alllabels if x[-1] in ['P','D']])

	
	return [posi,end, frontposi,clipsize]



def findbreak(t0, quality0,chrom0):

	alltags=[bin(int(x))[2:][::-1] for x in list(t0[1])]

	alltags=[str(x)+''.join(['0' for i in xrange(12-len(x))]) for x in alltags]

	quality=list(t0[4])

	cigar=list(t0[5])

	goodreads=[i for i,x in enumerate(alltags) if x[8]=='0' and x[-1]=='0' and x[0]=='1' and x[1]=='1' and x[2]=='0' and x[3]=='0'  and int(quality[i])>quality0 ]
	
	t0=t0.iloc[goodreads]
	
	cigar=list(t0[5])
	
	posis=list(t0[3])
	
	
	clip_position=[readcigar(x,int(posis[i]))+[int(quality[i])] for i,x in enumerate(cigar)]
	
	print 'findclip ',len(clip_position)
	
	out=pd.DataFrame.from_records(clip_position).to_csv('{:s}_findclip.csv'.format(chrom0),mode='a',header=None, index=False, sep=',')
	
	return clip_position
	
	

def main(inputfile):
	

	#allchroms=['chr5', 'chr6', 'chr12', 'chr16', 'chr7', 'chr10', 'chr14', 'chr15', 'chr8', 'chr19', 'chr9', 'chrY', 'chr18', 'chrM', 'chr17', 'chr3', 'chr13', 'chr20', 'chr11', 'chr2', 'chr1', 'chr4', 'chr21']

	allchroms=['chrX']
	
	try:
		os.mkdir('temp')	
	except:
		pass

	for chrom in allchroms:
		
		chromfile='temp/'+inputfile.split('/')[-1]+'_'+chrom
	
		if chrom != 'chr5':	
			os.system('samtools view {:s} {:s} > {:s}'.format(inputfile, chrom, chromfile))
	
		t = pd.read_csv(chromfile, sep='\t',dtype='str', chunksize=100000,names=range(30))
	
		for chunk in t:
		
			findbreak(chunk, 20, chrom)
	
	
if __name__=='__main__':
	
	inputfile=''
	opts,args=getopt.getopt(sys.argv[1:],"f:")
	for op, value in opts:
		if op=='-f':
			inputfile=value
	
	main(inputfile)



