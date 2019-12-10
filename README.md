# short-reads-clipping-detector

//  Created by Walfred MA in 2018, wangfei.ma@ucsf.edu.
//  Copyright © 2018 UCSF-Kwoklab. All rights reserved.

Finding positions of clippings from short reads alignments in bam files.

It has three pieces of scripts:

1. findclip1.py
2. findclip2.py
3. findclip3.py

findclip1.py:
It is used to find and record all positions of short reads clippings and their qualities.

Usage example: Python2.7 -f inputbamfile.bam

outputs to: inputbamfile_chr*_findclips.csv


findclip2.py:
It is used to find clipping rates on all reference positions.

Usage example: Python2.7 -f inputbamfile_chr*_findclips.csv

outputs: inputbamfile_chr*_findclips2.csv


findclip3.py:
It is used find all positions with high clippings rate on designated regions.

Usage example: Python2.7 -f SV.files -r inputbamfile_chr*_findclips.csv

outputs: SV.files_withclipping
