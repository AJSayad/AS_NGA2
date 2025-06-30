#clean files
echo "BASH - Cleaning old files."
rm -rf ensight/ monitor/
rm -rf singlephase_ensight/ singlephase_monitor/
rm -rf singlephase_log.out multiphase_log.out
rm -rf *.dat
rm -rf log.out

#reset extraction flag so that singlephase sim will run
echo "BASH - Reset extraction flag."
sed -i 's/Profile extraction flag: false/Profile extraction flag: true/' input

#reset multiphase timestep size so that it can be replaced to match final timestep of singlephase sim
echo "BASH - Reset multiphase max timestep size."
sed -i '/Multiphase max timestep size:/ s/:.*$/: @replace@/' input

#source config to set envoirnment variables
echo "BASH - Sourcing config.sh file."
source config.sh

#run singlephase simulation
echo "BASH - Running singlephase simulation in serial."

#update processors for singlephase serial simulation
sed -i '/Partition :/ s/:.*$/: 1 1 1/' input
mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i ${file_in} -v 2 > ${singlephase_file_out} 

#update extraction flag so that the multiphase sim runs
sed -i 's/Profile extraction flag: true/Profile extraction flag: false/' input

#update timestep in input file
cd monitor/
export filename='simulation'
export file_col=3

tail -1 simulation | while read line; do
    value=$(echo $line | cut -d " " -f $file_col)
    echo "BASH - Max timestep size for multiphase simualtion" $value

    sed -i "s/@replace@/$value/g" ../input 
done

cd ..

#update processors for multiphase based on config file
echo 'BASH - Updating domain decomposition for multiphase.'
sed -i "s/Partition : 1 1 1/Partition : $nprocx $nprocy $nprocz/g" input 

#store singlphase files
echo "BASH - Moving singlephase monitor and ensight directories."
mv ensight/ singlephase_ensight/
mv monitor/ singlephase_monitor/

#run multiphase simulation
echo "BASH - Running multiphase simulation."
mpiexec -n ${ntasks} ./nga.dp.gnu.opt.mpi.exe -i ${file_in} -v 2 > ${multiphase_file_out}
