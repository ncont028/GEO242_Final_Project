#!/bin/bash

data_dir='Hayward_RE_catalogs/'
file_str='RE_catalog_'
file_ext='.txt'

for g in 001 002 003 004 005 006 007 008 009 010 011 012 013 014
do
	gfile=$data_dir$file_str$g$file_ext
	echo $gfile
	#grep -vn : RE_catalog_001.txt | awk -F: '{print $1}'>>g1_empty.csv
	outfile='g'$g'_empty.csv'
	#echo $outfile
	grep -vn : $gfile | awk -F: '{print $1}' >> $outfile
done
