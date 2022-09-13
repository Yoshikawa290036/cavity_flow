#!/usr/bin/bash
ffmpeg -framerate 20 -start_number 0 -i image/xyuvp%03d000.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vframes 1001 -vcodec libx264 -pix_fmt yuv420p -r 60 cavty_flow.mp4 -y
