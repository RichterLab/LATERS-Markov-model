# LATERS-Markov-model
LAgrangian Transport Eulerian Reaction Spatial (LATERS) Markov model for reactive transport

MATLAB scripts for the LATERS Markov Model
_________________________________________________________________________________________________

1. Fully Resolved model - Run master1.m to run the fully resolved model.

master1.m

-Specifies model file name.

-Runs rxnrw_fluxweight.m.

rxnrw_fluxweight.m

-Runs the fully resolved model. 

-Requires the flow field file, vel_field1.mat.

-Runs rk1.m, reflect.m, and rxn_cycle.m.

rk1.m

-Moves particles by random walk.

-Runs u_vel.m and v_vel.m to interpolate particle velocities from flow field.

u_vel.m

-Interpolates particle velocities in x direction from flow field.

v_vel.m

-Interpolates particle velocities in x direction from flow field.

reflect.m

-Reflects particles off top and bottom boundaries if a particle exits either boundary. 

rxn_cycle.m

-Identifies which particles (if any) react with one another, removes A and B particles that react, and creates C particles at reaction locations.

_________________________________________________________________________________________________

2. Parameterize LATERS Markov model – Run master2.m to run LATERS Markov model parameterization then run getTM.m to obtain the transition matrix needed to run the LATERS Markov model.

master2.m

-Specifies file name.

-Runs parameterize_fluxweight.m.

parameterize_fluxweight.m

-Tracks particles as they travel across two SMM cells. Records travel times and other necessary information to run the LATERS Markov model.

-Requires the flow field file, vel_field1.mat. 

-Runs rk1.m.

rk1.m

-Moves particles by random walk.

-Runs u_vel.m and v_vel.m to interpolate particle velocities from grid.

u_vel.m

-Interpolates particle velocities in x direction from flow field.

v_vel.m

-Interpolates particle velocities in x direction from flow field.

getTM.m

-uses parameterize_fluxweight.m output (parameterAB_LATERS_final.mat) to obtain the transition matrix needed to run the LATERS model. Saves output as TMdata.mat.

_________________________________________________________________________________________________

3. LATERS Markov model - Run master1.m to run the LATERS Markov model.

master1.m

-Specifies model file name.

-Runs run_LATERSAB_Method1.m.

-Requires the transition matrix file TMdata.mat. This file contains contains all information needed to run the LATERS Markov model including the transition matrix, number of particles, cell length, dispersion coefficient, domain size, velocity field, etc. Created with getTM.m in the model parameterization step.

run_LATERSAB_Method1.m

-Runs the LATERS Markov model by implementing transport with a Spatial Markov model and Eulerian reaction.
_________________________________________________________________________________________________

4. LATERS model where nearby A/B particle pairs are selected for reaction - Run master1.m to run the LATERS Markov model with nearby particle killing.

master1.m

-Specifies model file name.

-Runs run_LATERSAB_Method1_pkillnearby.m

-Requires the transition matrix file TMdata.mat. This file contains contains all information needed to run the LATERS Markov model including the transition matrix, 
number of particles, cell length, dispersion coefficient, domain size, velocity field, etc. Created with getTM.m in the model parameterization step.

run_LATERSAB_Method1_pkillnearby.m

-This is the same as the original LATERS Markov model, but A and B particles are selected for reaction based on how close they are to each other within the reaction cell. 

_________________________________________________________________________________________________

Additional files:

perm_orig1.mat = permeability field

vel_field1.mat = velocity field

TMdata.mat = transition matrix

