#!/bin/bash

# to extract files use tar --zstd -xvf 'archive.tar.zst' 

# compress smesh folder
tar --zstd -cvf droplet_plic.tar.zst droplet_plic

# move to field data folder
cd ShockDrop/

tar --zstd -cvf RHOG.tar.zst RHOG/
tar --zstd -cvf RHOL.tar.zst RHOL/
tar --zstd -cvf VOF.tar.zst VOF/
tar --zstd -cvf label.tar.zst label/
tar --zstd -cvf PG.tar.zst PG/
tar --zstd -cvf PL.tar.zst PL/
tar --zstd -cvf IG.tar.zst IG/
tar --zstd -cvf IL.tar.zst IL/
tar --zstd -cvf TG.tar.zst TG/
tar --zstd -cvf TL.tar.zst TL/
tar --zstd -cvf Mach.tar.zst Mach/
tar --zstd -cvf beta.tar.zst beta/
tar --zstd -cvf visc.tar.zst visc/
tar --zstd -cvf sound_speed.zst sound_speed/
tar --zstd -cvf velocity.tar.zst velocity/
