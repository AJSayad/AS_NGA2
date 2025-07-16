ffmpeg -framerate 12 -i schlieren_%03d.png -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p schlieren.mp4
