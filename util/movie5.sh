#!/bin/bash
#$ -cwd
#$ -o movstdout.txt
#$ -e movstderr.txt
#$ -M sjh518@york.ac.uk
#$ -l h_rt=120:00:00
date

module load ffmpeg/2.2.4


for i in `seq 1 5`;
do
	ffmpeg -y -framerate 24 -i out$i/sppframe%04d000.png -vf \
	"scale=iw*16:ih*16" ${PWD##*/}_${i}spp.mp4
        
	ffmpeg -y -framerate 24 -i out$i/lenframe%04d000.png -vf \
	"scale=iw*16:ih*16" ${PWD##*/}_${i}len.mp4
done
