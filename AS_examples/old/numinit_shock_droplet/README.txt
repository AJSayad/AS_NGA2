Numerically initialized shock droplet simulation. This runs a 1D shocktube simluation, extracts the shock profile and saves it to files with .dat extensions, and then runs a 2D multiphase simulation where the numerical shock profile is used to initialize the shock. 

How to Run:
1. Set up input file and set how far the singlephase shock should travel before writing the profile to .dat files.
2. Set up config.sh file with desired number of processors for running in parallel.
3. execute in terminal chmod +x config.sh run_sim.sh
4. To run full simulation (with numerically initialized multiphase sim) exectute in terminal ./run_sim.sh
