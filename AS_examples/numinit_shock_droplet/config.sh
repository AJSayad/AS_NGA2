export nprocx=4 #number of processors in x, y, and z
export nprocy=1
export nprocz=1
export ntasks=$((nprocx*nprocy*nprocz))

echo "BASH - Number of processors x: " $nprocx
echo "BASH - Number of processors y: " $nprocy
echo "BASH - Number of processors z: " $nprocz

export file_in=input #set input file
echo "BASH - Name of input file: " $file_in

export singlephase_file_out=singlephase_log.out #set singlephase log output file
echo "BASH - Singlephase output log: " $singlephase_file_out

export multiphase_file_out=multiphase_log.out #set multiphase log output file
echo "BASH - Multiphase output log: " $multiphase_file_out
