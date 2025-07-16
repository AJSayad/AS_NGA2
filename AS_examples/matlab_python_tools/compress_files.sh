#!/bin/bash

# compress smesh folder
tar -cvf droplet_smesh.tar.gz droplet_smesh

# move to field data folder
cd ShockDroplet/

tar -cvf Bulkmod.tar.gz	Bulkmod/
tar -cvf Density.tar.gz Density/
tar -cvf GP.tar.gz GP/
tar -cvf Grho.tar.gz Grho/
tar -cvf GrhoE.tar.gz GrhoE/
tar -cvf LP.tar.gz LP/
tar -cvf Lrho.tar.gz Lrho/
tar -cvf LrhoE.tar.gz LrhoE/
tar -cvf P.tar.gz P/
tar -cvf Tmptr.tar.gz Tmptr/
tar -cvf velocity.tar.gz velocity/
tar -cvf VOF.tar.gz VOF/


