#!/bin/bash


#This generates the movie
#ffmpeg -framerate 50 -i lenframe%05d00.png -c:v libx264 -vf "transpose=0,scale=iw*16:ih*16" -r 30 -pix_fmt yuv420p lenout100.mp4


#ALTERNATIVES
#The above command doesn't work on YARCC because libx264 isn't available - alternatives are: 

#ffmpeg -framerate 24 -i lenframe%05d00.png -vf "transpose=0,scale=iw*16:ih*16" output.mp4


ffmpeg -framerate 24 -i lenframe%04d000.png -vf "transpose=0,scale=iw*16:ih*16" output.mp4


echo "YARCC note: If there's been any problems, use 'module load ffmpeg/2.2.4' "
