%rw_fvm_fluxweight.m
% 
% clear all
% clc
% close all

function rxnrw_fluxweight(filename)

global V D

load vel_field1.mat

szV = size(V.x); %size of flow field
normV = mean(mean(V.x)); %mean velocity in x direction of flow field as is
meanu = 1; %desired mean velocity in x direction
V.x = meanu*V.x/normV; %rescale u to have <u> = meanu
V.y = meanu*V.y/normV; %rescale v in the same way as u

N = 5e5; %Number of A/B Particles (2N = total number of particles)
Nnow = N; %number of particles at current time step (will change with reactions)
D = 0.01; %Diffusion coefficient
k = 10; %reaction rate coefficient
xinit = 3200; %initial x position
lambda = 2; %correlation length
Lcell = 24*lambda; %SMM cell length
ABwidth = 10*Lcell; %width of area covered by the A and B particles at t=0

Nsteps = 64e4; %total number of steps
dt = 0.01; %T/total;    %Timestep size
tmax = Nsteps*dt; %time at end of simulation
time = dt:dt:tmax; 

Ly = szV(2); %domain length in y direction
Lx = szV(1); %domain length in x direction
dx = 1; %grid cell size in x direction
dy = 1; %grid cell size in y direction
Agrid = dx*dy; %grid cell area

Ly_part = 200; %length of y over which particles are present at start
ymin = 0.45*Ly; %minimum y value at start
ymax = 0.55*Ly; %maximum y value at start

%C0 = 1; %total mass of all the particles
%A0 = ABwidth*Ly_part; %initial area covered by A/B particles
Mtot = 1500; %total mass in system
mp = Mtot/N; %mass of each particle


%SET UP FLUX WEIGHTED INITIAL CONDITION
Npts = 1e3; %number of points at which we want particles placed in initial condition
yvec = linspace(ymin,ymax,Npts); %y position vector to use for initial condition
xvec = linspace(xinit,xinit+2*ABwidth,Npts); %x position vector to use for initial condition

uvec = zeros(Npts,Npts);
vvec = zeros(Npts,Npts);

for yidx = 1:Npts

    for xidx=1:Npts

        %find out which gridbox each x,y location is in
        x1 = floor(xvec(xidx));
        y1 = floor(yvec(yidx));

        x2 = x1+1;
        y2 = y1+1;

        uvec(xidx,yidx) = V.x(x1,y1)+(xvec(xidx)-x1)./(x2-x1).*(V.x(x2,y1)-V.x(x1,y1));
        vvec(xidx,yidx) = V.y(x1,y1)+(yvec(yidx)-y1)./(y2-y1).*(V.y(x1,y2)-V.y(x1,y1));


    end
end


absv = sqrt((uvec.^2)+(vvec.^2)); %absolute velocity at inlet

%Flux weight A
absvA = absv(Npts/2+1:end,:);
sz_absvA = size(absvA);
sum_absvA = sum(sum(absvA)); %sum of velocity values along y at xinit
flux_weightedA = round(N*absvA/sum_absvA); %flux weighted initial condition (round to get integer # particles)
diff_fluxA = sum(sum(flux_weightedA))-N; %correct to make sure that we have the right # of particles
idx_correctA = randi(sz_absvA(1)*sz_absvA(2),1,abs(diff_fluxA)); %create random indices to place extra/remove extra particles
for aa = 1:abs(diff_fluxA) %need loop because duplicates likely exist in idx_correctA
    flux_weightedA(idx_correctA(aa)) = flux_weightedA(idx_correctA(aa))+1; 
end

%Flux weight B
absvB = absv(1:Npts/2,:);
sz_absvB = size(absvB);
sum_absvB = sum(sum(absvB)); %sum of velocity values along y at xinit
flux_weightedB = round(N*absvB/sum_absvB); %flux weighted initial condition (round to get integer # particles)
diff_fluxB = sum(sum(flux_weightedB))-N; %correct to make sure that we have the right # of particles
idx_correctB = randi(sz_absvB(1)*sz_absvB(2),1,abs(diff_fluxB)); %create random indices to place extra/remove extra particles
for bb = 1:abs(diff_fluxB) %need loop because duplicates likely exist in idx_correctA
    flux_weightedB(idx_correctB(bb)) = flux_weightedB(idx_correctB(bb))+1; 
end


xgridA = repmat(xvec(Npts/2+1:end),Npts,1)';
xgridB = repmat(xvec(1:Npts/2),Npts,1)';
ygrid = repmat(yvec',1,Npts/2)';


xA = zeros(1,N);
yA = zeros(1,N);
    idx1A = 1;
    idx2A = flux_weightedA(1,1);
    
xB = zeros(1,N);
yB = zeros(1,N);
    idx1B = 1;
    idx2B = flux_weightedB(1,1);

for cc = 1:sz_absvB(1)
    
    for dd = 1:sz_absvB(2)
    
        %Determine xA, yA, xB, and yB
        if flux_weightedA(cc,dd)>0
            xA(idx1A:idx2A) = xgridA(cc,dd);
            yA(idx1A:idx2A) = ygrid(cc,dd);
        end
        if flux_weightedB(cc,dd)>0
            xB(idx1B:idx2B) = xgridB(cc,dd);
            yB(idx1B:idx2B) = ygrid(cc,dd);
        end
    
        %update indices for flux weighting
        if dd<sz_absvB(2)
            
              if flux_weightedA(cc,dd+1)>0
                idx1A = idx2A+1;
                idx2A = idx2A+flux_weightedA(cc,dd+1);
              end
              
              if flux_weightedB(cc,dd+1)>0
                idx1B = idx2B+1;
                idx2B = idx2B+flux_weightedB(cc,dd+1);
              end 
        end
        
        if dd==sz_absvB(2) & cc<sz_absvB(1)
              if flux_weightedA(cc+1,1)>0
                idx1A = idx2A+1;
                idx2A = idx2A+flux_weightedA(cc+1,1);
              end
              
              if flux_weightedB(cc+1,1)>0
                idx1B = idx2B+1;
                idx2B = idx2B+flux_weightedB(cc+1,1);
              end    
        end

    end
            
end

%store vector of the initial particle positions
xAinit = xA;
yAinit = yA;
xBinit = xB;
yBinit = yB;

xC = []; %product C is initially empty, no reactions occurred before t = 0
yC = [];

%vectors to store particles 10 most recent locations before reaction
xArecent=[zeros(N,9),xA'];
xBrecent=[zeros(N,9),xB'];
yArecent=[zeros(N,9),yA'];
yBrecent=[zeros(N,9),yB'];

react_loc   = zeros(N,3);           %vector to store locations and timestep of reactions
countreact  = 0;                    %initially no particles have reacted
xArecent_loc = zeros(N,10);         %matrix that stores 10 most recent particle locations before reaction
yArecent_loc = zeros(N,10);         %matrix that stores 10 most recent particle locations before reaction
xBrecent_loc = zeros(N,10);         %matrix that stores 10 most recent particle locations before reaction
yBrecent_loc = zeros(N,10);         %matrix that stores 10 most recent particle locations before reaction


meanCA = zeros(1,Nsteps);
meanCB = zeros(1,Nsteps);
meanCC = zeros(1,Nsteps);


for kk=1:Nsteps
    
     %kk
     
     Nnow = length(xB); %current number of A and B particles
     NC = N-Nnow; %current number of C particles
     
     Pr = k*mp*dt; %probability of reaction given collocation (using max mass of B particle)
     P = 0.000001; %probability of reaction
     r_react = sqrt(-8*D*dt*log(8*pi*D*dt*P/Pr)); %radius for particle search in reaction step
     
     idxout = find([xB,xA]>Lx);
     Nout = length(idxout); %number of particles that have crossed the right boundary (shouldn't happen)

     if Nout>0 %stop if any particles have crossed right boundary of flow field
         break
     end
     
     %move by random walk
     [xA, yA] = rk1 (xA, yA, dt); 
     [xB, yB] = rk1 (xB, yB, dt); 
     [xC, yC] = rk1 (xC, yC, dt); 
     
     
     %reflect particles if they exit boundary (shouldn't be a problem)
     reflect
        
        %reset count grid for concentration grid calculation
        mAcount = zeros(szV);
        mBcount = zeros(szV);
        mCcount = zeros(szV);
        %bilinear interpolation to grid for mixing metric calculation
        for mm=1:Nnow

            %PREPARE CONCENTRATION GRID
            %identify which mixing grid box each particle is in
            rnddownxA = floor(xA(mm)/dx); %don't use xtemp here, we want to use actual x value (not mod(x,Lx))
            rnddownyA = floor(yA(mm)/dy);
            rnddownxB = floor(xB(mm)/dx); %don't use xtemp here, we want to use actual x value (not mod(x,Lx))
            rnddownyB = floor(yB(mm)/dy);
            
            mAcount(rnddownxA,rnddownyA) = mAcount(rnddownxA,rnddownyA)+1;
            mBcount(rnddownxB,rnddownyB) = mBcount(rnddownxB,rnddownyB)+1; %add one to count grid

        end
        
        for nn=1:NC

            %PREPARE CONCENTRATION GRID
            %find out which mixing grid box each particle is in
            rnddownxC = floor(xC(nn)/dx); %don't use xtemp here, we want to use actual x value (not mod(x,Lx))
            rnddownyC = floor(yC(nn)/dy);
            
            mCcount(rnddownxC,rnddownyC) = mCcount(rnddownxC,rnddownyC)+1;

        end
        
        CA = mAcount*mp/Agrid; %concentration of A grid
        CB = mBcount*mp/Agrid; %concentration of B grid
        CC = mCcount*2*mp/Agrid; %concentration of C grid
        meanCA(kk) = mean(mean(CA));
        meanCB(kk) = mean(mean(CB));
        meanCC(kk) = mean(mean(CC));
        
        posA = [xA;yA]; %gives position matrix of x in row 1, y in row 2
        posB = [xB;yB];

        xyA = transpose(posA); %first column x, second column y
        xyB = transpose(posB);

        %stop simulation if all particles have reacted
        if isempty(xyA)==1
            break
        end


        MdlKDT = KDTreeSearcher(xyB);   %create kd tree of B particles

        clear IdxTot
        
        IdxTot = rangesearch(MdlKDT,xyA,r_react); %find indices of B particles within a distance r of the A particles
        
        rxn_cycle
     
     
%      figure(1)
%      plot(xB,yB,'r.')
%      hold on
%      plot(xA,yA,'b.')
%      plot(xC,yC,'g.')


        if mod(kk,500)==0
            
            meanCAnow = meanCA(kk);
            meanCBnow = meanCB(kk);
            meanCCnow = meanCC(kk);

            save([filename,'_pos_time',int2str(kk),'.mat'],'xA','yA','xB',...
                 'yB','xC','yC','meanCAnow','meanCBnow','meanCCnow');

        end
        

end

save([filename,'.mat'])

end
