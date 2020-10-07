%LATERS Markov Model - use Method 1 to downscale particle locations
%instead of removing particle mass on each particle evenly across a cell,
%remove particles based on their probability of reaction. 

% clear
% clc
% close all

function run_LATERSAB_Method1(filename)

load('TMdata.mat')
%load('vel_field.mat')

%N = 1e6; %total number of particles
N_AB = N/2; %number of A/B particles (N/2 since half are A and half are B particles)
%lambda = 8; %correlation length of log permeability field
%Lcell = 6*lambda; %length of cell
%meanu = 1;
%tmin = lambda/meanu; %mininum time at which we can expect upscaled model to be accurate

k = 10; %reaction rate

dt = 0.1; %time step at which we calculate reactions
t = dt:dt:6400; %times at which to measure reactions 

Nsteps = length(t); %number of steps at which we calculate reactions

rng('shuffle')

%%% GRID INFO %%%
Ly = szV(2); %Length of domain in y direction
Lx = szV(1)-1; %Length of domain in x direction
dx = 2; %grid cell size in x direction for concentration grid (Lx/dx must be an integer)
dy = 2; %grid cell size in y direction for concentration grid (Ly/dy must be an integer)
Agrid = dx*dy; %grid cell area
% Ngridx = Lx/dx; %number of grid cells in x direction
% Ngridy = Ly/dy; %number of grid cells in y direction
% xgrid = 0:dx:Lx;
% ygrid = 0:dy:Ly;

rxncell_x = 48; %x length of volume averaged rxn cell
rxncell_y = 48; %y length of volume averaged rxn cell
Arxncell = rxncell_x*rxncell_y; %volume of rxn cell

Ly_part = max(y0)-min(y0); %length of y over which particles are present at start
ymin = min(y0); %minimum y value at start

Mtot = 1500; %total mass of A (same as total mass of B)
mp = Mtot/N_AB; %mass of each particle

xp = x0+Lcell; %particle x locations after one SMM step
yp = y1_inlet; %particle y locations after one SMM step

tp = tau1; %vector of particle times after one SMM step, move with tau1 in first step
tau = tau1; %vector of particle taus (will be updated each SMM step)
bin_now = bins_tau1'; %vector of particle bin numbers for this SMM step

TM_prob = [zeros(numbins,1),cumsum(TM,2)]; %add zeros for new bin search (represents left edge of bin probabilities)

bin_prev = zeros(1,N); %vector of bin numbers for all particles in previous SMM step
xp_prev = x0; %particle x locations from previous SMM step
yp_prev = y0_inlet; %particle y locations at outlet of previous cell (inlet of current cell)
tp_prev = zeros(1,N); %vector of total particle times to cross previous cell (left edge of new cell)
tau_prev = zeros(1,N); %vector of particle travel times through previous cell
yp0 = y0_inlet; %particle's y value at left edge of previous cell
yp1 = y1_inlet; %particle's y value at right edge of previous cell
dyp = y1_inlet-y0_inlet; %change in y position across first SMM cell

meanCB_tot = zeros(1,Nsteps); %vector of CB averaged over entire domain at each reaction step
meanCA_tot = zeros(1,Nsteps); %vector of CD averaged over entire domain at each reaction step (should be zero)
mB_tot = zeros(1,Nsteps); %sum of all the mass of B in the system at each reaction step
mA_tot = zeros(1,Nsteps); %sum of all the mass of D in the system at each reaction step

mpB = [mp*ones(1,N_AB),zeros(1,N_AB)]; %mass of each particle
mpA = [zeros(1,N_AB),mp*ones(1,N_AB)]; %mass of each particle

for kk = 1:Nsteps
    
    %kk
    
    tnow = t(kk);
    
    %NOTE: t(kk) falls within the bin labeled "bin_now" at each step. We
    %store the information from the previous step because it is needed for
    %the downscaling calculations. The amount of time spent by the particle
    %in the current cell is t(kk)-tp_prev.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% TRANSPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for pp = 1:N %loop through all particles
        
        while tp(pp)<tnow %while the particle's time is less than the time 
                          %at which we are approximating mixing/reaction
            
            %store tau, x position, and total particle time from previous step              
            tau_prev(pp) = tau(pp);
            xp_prev(pp) = xp(pp);
            tp_prev(pp) = tp(pp); 
            bin_prev(pp) = bin_now(pp);
    
            %%Determine new bin for particle in this step
            xi = rand; %random number from U(0,1) to determine new bin number
    
            for ii = 1:numbins

                if xi<TM_prob(bin_prev(pp),ii+1) && xi>=TM_prob(bin_prev(pp),ii)
          
                    bin_now(pp) = ii;

                        %%Draw a random time from travel times in the bin to be each particle's tau value
                        if ii == 1
                            %Draw random tau from travel times tau1 that were in state 1 (times_1)
                            tau(pp) = times_1(randi([1,ntimes1],length(pp),1));
                        end

                        if ii == 2
                            tau(pp) = times_2(randi([1,ntimes2],length(pp),1)); 
                        end

                        if ii == 3
                            tau(pp) = times_3(randi([1,ntimes3],length(pp),1)); 
                        end

                        if ii == 4
                            tau(pp) = times_4(randi([1,ntimes4],length(pp),1));
                        end

                        if ii == 5
                            tau(pp) = times_5(randi([1,ntimes5],length(pp),1));
                        end

                        if ii == 6
                            tau(pp) = times_6(randi([1,ntimes6],length(pp),1));
                        end

                        if ii == 7
                            tau(pp) = times_7(randi([1,ntimes7],length(pp),1)); 
                        end

                        if ii == 8
                            tau(pp) = times_8(randi([1,ntimes8],length(pp),1));
                        end

                        if ii == 9
                            tau(pp) = times_9(randi([1,ntimes9],length(pp),1));
                        end

                        if ii == 10
                            tau(pp) = times_10(randi([1,ntimes10],length(pp),1));
                        end

                        if ii == 11
                            tau(pp) = times_11(randi([1,ntimes11],length(pp),1));
                        end

                        if ii == 12
                            tau(pp) = times_12(randi([1,ntimes12],length(pp),1));
                        end

                        if ii == 13
                            tau(pp) = times_13(randi([1,ntimes13],length(pp),1));  
                        end

                        if ii == 14
                            tau(pp) = times_14(randi([1,ntimes14],length(pp),1)); 
                        end

                        if ii == 15
                            tau(pp) = times_15(randi([1,ntimes15],length(pp),1));   
                        end

                        if ii == 16
                            tau(pp) = times_16(randi([1,ntimes16],length(pp),1));
                        end

                        if ii == 17
                            tau(pp) = times_17(randi([1,ntimes17],length(pp),1)); 
                        end

                        if ii == 18
                            tau(pp) = times_18(randi([1,ntimes18],length(pp),1));
                        end

                        if ii == 19
                            tau(pp) = times_19(randi([1,ntimes19],length(pp),1));
                        end

                        if ii == 20
                            tau(pp) = times_20(randi([1,ntimes20],length(pp),1));
                        end
                end

            end

            %move particles forward in space and time
            xp(pp) = xp_prev(pp)+Lcell; 
            tp(pp) = tp_prev(pp)+tau(pp);
 
        end
        
    end
    
    %downscaled x and y values for mixing metric calculation
    x_ds = xp_prev+(tnow-tp_prev)./tau*Lcell;
    y_ds = y0_inlet; %Method 1: each particle maintains its initial y position at all times
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% REACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    minx = min(x_ds);
    maxx = max(x_ds);
    miny = min(y_ds);
    maxy = max(y_ds);
    
    Nrxncells_x = ceil((maxx-minx)/rxncell_x); %number of reaction cells within grid in x direction
    Nrxncells_y = ceil((maxy-miny)/rxncell_y); %number of reaction cells within grid in y direction
    
    lx_grid = Nrxncells_x*rxncell_x; %length of grid in x direction (need to have integer number of cells)
    ly_grid = Nrxncells_y*rxncell_y; %length of grid in y direction (need to have integer number of cells)
    
    Ngridx = lx_grid/dx+1; %number of grid cells in x direction
    Ngridy = ly_grid/dy+1; %number of grid cells in y direction
    
    xgrid = minx:dx:minx+Ngridx*dx; 
    ygrid = miny:dy:miny+Ngridy*dy;
    
    
    mAgrid = zeros(Ngridx,Ngridy); %grid to calculate the mass of D
    mBgrid = zeros(Ngridx,Ngridy); %grid to calculate the mass of B
    %bilinear interpolation to grid for reaction calculation
    for mm=1:N %loop over all particles

        %PREPARE CONCENTRATION GRID
        %find out which grid box each particle is in
        rnddownx = floor((x_ds(mm)-minx)/dx)+1; %add one because index starts at 1 not 0
        rnddowny = floor((y_ds(mm)-miny)/dy)+1;

        mAgrid(rnddownx,rnddowny) = mAgrid(rnddownx,rnddowny)+mpA(mm); %add particle D mass to grid
        mBgrid(rnddownx,rnddowny) = mBgrid(rnddownx,rnddowny)+mpB(mm); %add particle B mass to grid

    end
    
    %DIVIDE THE DOMAIN INTO CELLS OF WIDTH Lcell TO CALCULATE THE AMOUNT OF REACTION
    %Record the following information in each of these reaction cells:
    meanCArecord_cell = zeros(Nrxncells_x,Nrxncells_y); %<C_A>, where <.> is an average over the reaction cell
    meanCBrecord_cell = zeros(Nrxncells_x,Nrxncells_y); %<C_B>, where <.> is an average over the reaction cell
    term1 = zeros(Nrxncells_x,Nrxncells_y); %term 1 in the reaction equation in each reaction cell
    term2 = zeros(Nrxncells_x,Nrxncells_y); %term 2 in the reaction equation in each reaction cell
    r = zeros(Nrxncells_x,Nrxncells_y); %r = k*(term1+term2) 

    %loop through every reaction cell and calculate the amount of reaction in each cell
    for cellnum_x = 1:Nrxncells_x
        for cellnum_y = 1:Nrxncells_y
        
            idx_min_x = 1+(rxncell_x/dx)*(cellnum_x-1); %left edge of current cell
            idx_max_x = idx_min_x+(rxncell_x/dx); %right edge of current cell
            idx_min_y = 1+(rxncell_y/dy)*(cellnum_y-1); %bottom edge of current cell
            idx_max_y = idx_min_y+(rxncell_y/dy); %top edge of current cell

            CB_cell = mBgrid(idx_min_x:idx_max_x-1,idx_min_y:idx_max_y-1)/Agrid; %grid of reactive solute concentration CB
            CA_cell = mAgrid(idx_min_x:idx_max_x-1,idx_min_y:idx_max_y-1)/Agrid; %grid of reactive solute concentration CA
            meanCA_cell = mean(mean(CA_cell));  %<C_A>
            meanCB_cell = mean(mean(CB_cell));  %<C_B>

            meanCArecord_cell(cellnum_x,cellnum_y) = meanCA_cell; %<C_A>, where <.> is an average over the reaction cell
            meanCBrecord_cell(cellnum_x,cellnum_y) = meanCB_cell; %<C_B>, where <.> is an average over the reaction cell

            CAf_cell = CA_cell-meanCA_cell; %C_A' in cell
            CBf_cell = CB_cell-meanCB_cell; %C_B' in cell
            term1(cellnum_x,cellnum_y) = meanCB_cell*meanCA_cell; %<C_B><C_A>
            term2_cell = CBf_cell.*CAf_cell; %C_B'C_A'
            term2(cellnum_x,cellnum_y) = mean(mean(term2_cell)); %<C_B'C_A'>

            r(cellnum_x,cellnum_y) = k*(term1(cellnum_x,cellnum_y)+term2(cellnum_x,cellnum_y)); %calculate reaction term

            dCB = -r(cellnum_x,cellnum_y)*dt; %total change in C_B in the reaction cell due to reaction

            if dCB>0
                dCB = 0; %we cannot gain mass in our system
            end

            idx_incell = find(x_ds>=xgrid(idx_min_x) & x_ds<xgrid(idx_max_x) & y_ds>=ygrid(idx_min_y) & y_ds<ygrid(idx_max_y)); %indices of particles in current cell
            idxBpart = find(mpB(idx_incell)>0); %determine which particles are B particles with mass remaining
            idxApart = find(mpA(idx_incell)>0); %determine which particles are A particles with mass remaining
            idxB_rxn = idx_incell(idxBpart); %idx of B particles available for rxn
            idxA_rxn = idx_incell(idxApart); %idx of A particles available for rxn
            NB_rxn = length(idxB_rxn); %number of B particles in cell
            NA_rxn = length(idxA_rxn); %number of B particles in cell


            %adjust particle B mass
            mBtot_incell = sum(mpB(idxB_rxn)); %total mass of B in rxn cell - same as sum(sum(mBcount(idx_min:idx_max,:)))
            mBnew_incell = mBtot_incell+dCB*Arxncell; %new mass of B in rxn cell
            diff_mB = mBtot_incell-mBnew_incell;
            
            mAtot_incell = sum(mpA(idxA_rxn)); %total mass of A in rxn cell - same as sum(sum(mAcount(idx_min:idx_max,:)))
            mAnew_incell = mAtot_incell+dCB*Arxncell; %new mass of A in rxn cell
            diff_mA = mAtot_incell-mAnew_incell;

            %P_B_rxn and P_A_rxn are not the same because there are
            %(likely) different numbers of A and B particles in each
            %reaction cell
            P_B_rxn = diff_mB/mBtot_incell; %probability of reaction for B particles
            P_A_rxn = diff_mA/mAtot_incell; %probability of reaction for A particles

            if P_B_rxn>0
                P_B_part = rand(NB_rxn,1);
                idxB_remove = find(P_B_part<P_B_rxn);
                P_A_part = rand(NA_rxn,1);
                idxA_remove = find(P_A_part<P_A_rxn);
                
                N_remove = min([length(idxA_remove),length(idxB_remove)]); %make sure we remove the same number of A and B particles
                
                mpB(idxB_rxn(idxB_remove(1:N_remove))) = 0;
                mpA(idxA_rxn(idxA_remove(1:N_remove))) = 0;
            end
        
        end
    end
    
    meanCB_tot(kk) = mean(mean(mBgrid/Agrid)); %mean concentration averaged over whole domain before rxn
    meanCA_tot(kk) = mean(mean(mAgrid/Agrid)); %mean concentration averaged over whole domain before rxn

    mB_tot(kk) = sum(mpB); %total mass of B in system
    mA_tot(kk) = sum(mpA); %total mass of D in system

    if kk<100 && mod(kk,5)==0
        save(['paper_latersAB_SMM_Method1_pkill_vcell_newgrid_dxdy8_rxndxdy48_pos_time',int2str(kk),'.mat'],'term1','term2','r','mpB','mpA','x_ds','y_ds','meanCArecord_cell','meanCBrecord_cell')
    end
    
    if mod(kk,100)==0     
        save(['paper_latersAB_SMM_Method1_pkill_vcell_newgrid_dxdy8_rxndxdy48_pos_time',int2str(kk),'.mat'],'term1','term2','r','mpB','mpA','x_ds','y_ds','meanCArecord_cell','meanCBrecord_cell')  
    end
  
end

save('paper_latersAB_SMM_Method1_pkill_vcell_newgrid_all.mat','mB_tot','mA_tot','meanCB_tot','meanCA_tot')

end
