#!/bin/bash
for file in *.fastq.gz
do
 echo "Unziping"
 gzip -d  "$file"
done

#Rename the files
for file in *.fastq
do
 echo "Renaming"
 newname=$(echo $file | cut -d_ -f1,2,5).fastq
 echo "Renaming $file as $newname"
# mv $file $newname 
done
##Script from HPC workshop 2 3/16/2023
##Updated on 6/7/2023 with troubleshooting with ChatGPT. 
#Various itterations offered were not exactly what
# I needed.-BAR

