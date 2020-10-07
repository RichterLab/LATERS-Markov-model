%parameterize
%RECORD TRAJECTORIES OF PARTICLES THROUGH CELL

function parameterize_fluxweight(filename)

global V D

load vel_field1.mat

szV = size(V.x);
normV = mean(mean(V.x));
meanu = 1;
V.x = meanu*V.x/normV; %rescale u to have meanu as mean velocity in x direction
V.y = meanu*V.y/normV; %rescale v in the same way as u

N = 5e5;      %Number of particles
D = 0.01;     %Constant dispersion coefficient
xinit = 3200; %Initial x position
lambda = 2;   %Correlation length
Lcell = 24*lambda; %Length of SMM cell
ABwidth = 10*Lcell; %Width of initial particle distribution

total = 10000000; %total number of steps - set very large to allow most 
                  %(if not all) particles to travel two cell lengths
dt = 0.001; %Timestep size
%T = total*dt;
%time = dt:dt:T;

Lx = szV(1); %domain length in x direction
Ly = szV(2); %domain length in y direction

ymin = 0.45*Ly; %minimum y value at t=0
ymax = 0.55*Ly; %maximum y value at t=0
Ly_part = ymax-ymin; %length of y over which particles are present at start

%C0 = 1; %total mass of all the particles
%A0 = ABwidth*Ly_part; %initial area covered by A/B particles
Mtot = 1500; %total mass in the system
mp = Mtot/N; %mass of each particle

%SET UP FLUX WEIGHTED INITIAL CONDITION
Npts = 1e3; %number of points to place particles in initial positions
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
for bb = 1:abs(diff_fluxB) %need loop because duplicates likely exist in idx_correctB
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

%Complete parameterization for both A and B particles
x = [xB,xA]; %all particle (both A and B) y positions
y = [yB,yA]; %all particle (both A and B) x positions
x0 = x; %initial x positions
y0 = y; %initial y positions
N = length(x); %total number of particles (both A and B)

times1 = zeros(1,N); %vector to track travel times through first SMM cell
times2 = zeros(1,N); %vector to track travel times through first SMM cell
y0_inlet = y0;
y1_inlet = zeros(1,N);

%Keep track of particle x positions at 20 locations across the SMM cell
traj1 = zeros(1,N);
traj2 = zeros(1,N);
traj3 = zeros(1,N);
traj4 = zeros(1,N);
traj5 = zeros(1,N);
traj6 = zeros(1,N);
traj7 = zeros(1,N);
traj8 = zeros(1,N);
traj9 = zeros(1,N);
traj10 = zeros(1,N);
traj11 = zeros(1,N);
traj12 = zeros(1,N);
traj13 = zeros(1,N);
traj14 = zeros(1,N);
traj15 = zeros(1,N);
traj16 = zeros(1,N);
traj17 = zeros(1,N);
traj18 = zeros(1,N);
traj19 = zeros(1,N);

%Keep track of particle y positions at 20 locations across the SMM cell
y_traj1 = zeros(1,N);
y_traj2 = zeros(1,N);
y_traj3 = zeros(1,N);
y_traj4 = zeros(1,N);
y_traj5 = zeros(1,N);
y_traj6 = zeros(1,N);
y_traj7 = zeros(1,N);
y_traj8 = zeros(1,N);
y_traj9 = zeros(1,N);
y_traj10 = zeros(1,N);
y_traj11 = zeros(1,N);
y_traj12 = zeros(1,N);
y_traj13 = zeros(1,N);
y_traj14 = zeros(1,N);
y_traj15 = zeros(1,N);
y_traj16 = zeros(1,N);
y_traj17 = zeros(1,N);
y_traj18 = zeros(1,N);
y_traj19 = zeros(1,N);
    

for idx=1:total
    
     %idx
     
     idxmove = find(x<(x0+2*Lcell+2)); %idx of particles that have not yet 
                                       %traveled two SMM cell lengths - we
                                       %stop moving particles after they
                                       %travel 2Lcell from x0
     Nmove = length(idxmove); %number of particles that continue to move

     if Nmove==0 %stop if all particles have passed last BTC location (100L)
         break
     end
     
     xtemp = x(idxmove);
     ytemp = y(idxmove);
        
     [xnewtemp, ynewtemp] = rk1 (xtemp, ytemp, dt);
     
     xnew = x;
     ynew = y;
     xnew(idxmove) = xnewtemp;
     ynew(idxmove) = ynewtemp;
    
        %reflect particles off top and bottom boundary
            %top boundary
            out=find(y>Ly);    
            y(out)=2*Ly-y(out); 
            clear out 

            %bottom boundary
            out=find(y<1);    
            y(out)=2-y(out);
            clear out

     
     %measure particle trajectories at 20 locations across SMM cell
     idxbtc = find(xnew>(x0+0.05*Lcell) & x<(x0+0.05*Lcell));
     traj1(idxbtc) = idx*dt;
     y_traj1(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.1*Lcell) & x<(x0+0.1*Lcell));
     traj2(idxbtc) = idx*dt;
     y_traj2(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.15*Lcell) & x<(x0+0.15*Lcell));
     traj3(idxbtc) = idx*dt;
     y_traj3(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.2*Lcell) & x<(x0+0.2*Lcell));
     traj4(idxbtc) = idx*dt;
     y_traj4(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.25*Lcell) & x<(x0+0.25*Lcell));
     traj5(idxbtc) = idx*dt;
     y_traj5(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.3*Lcell) & x<(x0+0.3*Lcell));
     traj6(idxbtc) = idx*dt;
     y_traj6(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.35*Lcell) & x<(x0+0.35*Lcell));
     traj7(idxbtc) = idx*dt;
     y_traj7(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.4*Lcell) & x<(x0+0.4*Lcell));
     traj8(idxbtc) = idx*dt;
     y_traj8(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.45*Lcell) & x<(x0+0.45*Lcell));
     traj9(idxbtc) = idx*dt;
     y_traj9(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.5*Lcell) & x<(x0+0.5*Lcell));
     traj10(idxbtc) = idx*dt;
     y_traj10(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.55*Lcell) & x<(x0+0.55*Lcell));
     traj11(idxbtc) = idx*dt;
     y_traj11(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.6*Lcell) & x<(x0+0.6*Lcell));
     traj12(idxbtc) = idx*dt;
     y_traj12(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.65*Lcell) & x<(x0+0.65*Lcell));
     traj13(idxbtc) = idx*dt;
     y_traj13(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.7*Lcell) & x<(x0+0.7*Lcell));
     traj14(idxbtc) = idx*dt;
     y_traj14(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.75*Lcell) & x<(x0+0.75*Lcell));
     traj15(idxbtc) = idx*dt;
     y_traj15(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.8*Lcell) & x<(x0+0.8*Lcell));
     traj16(idxbtc) = idx*dt;
     y_traj16(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.85*Lcell) & x<(x0+0.85*Lcell));
     traj17(idxbtc) = idx*dt;
     y_traj17(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.90*Lcell) & x<(x0+0.90*Lcell));
     traj18(idxbtc) = idx*dt;
     y_traj18(idxbtc) = ynew(idxbtc);
     
     idxbtc = find(xnew>(x0+0.95*Lcell) & x<(x0+0.95*Lcell));
     traj19(idxbtc) = idx*dt;
     y_traj19(idxbtc) = ynew(idxbtc);
     
     
     %record travel times at breakthrough curve locations
     idxbtc = find(xnew>(x0+Lcell) & x<(x0+Lcell));
     times1(idxbtc) = idx*dt; %this is also traj20
     y1_inlet(idxbtc) = ynew(idxbtc); 

     idxbtc = find(xnew>(x0+2*Lcell) & x<(x0+2*Lcell));
     times2(idxbtc) = idx*dt;

            
     x = xnew;
     y = ynew;
    
end

save([filename,'_final.mat'])

end
