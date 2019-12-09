# short-reads-clipping-detector
Finding the position of clippings from short reads alignments in the bam files

There are three pieces of scripts:

1. findclip1.py
2. findclip2.py
3. findclip3.py

findclip1.py:
It is used to find and record all short reads clips positive and their qualities.

Usage example: Python2.7 -f inputbamfile.bam

outputs to: inputbamfile_chr*_findclips.csv


findclip2.py:
It is used to find clipping rates on all reference positions.

Usage example: Python2.7 -f inputbamfile_chr*_findclips.csv

outputs: inputbamfile_chr*_findclips2.csv


findclip3.py:
It is used to overlap clips on all structural variation from 10X table file.

Usage example: Python2.7 -f SV.files -r inputbamfile_chr*_findclips.csv

outputs: SV.files_withclipping
