#!/bin/bash


#This generates the movie
ffmpeg -framerate 50 -i lenframe%05d00.png -c:v libx264 -vf "transpose=0,scale=iw*16:ih*16" -r 30 -pix_fmt yuv420p lenout100.mp4

