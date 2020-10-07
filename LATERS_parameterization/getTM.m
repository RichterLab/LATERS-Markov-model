%Nikki's suggestion: create transition matrix with equiprobable bins for
%tau1 and use the same bins for tau2 (which will not be perfectly
%equiprobable). Sample times for SMM from tau1 since it is equiprobable.

numbins = 20;
N = 100000;

load('parameterAB_LATERS_final.mat')

tau1 = times1; %time to travel through cell 1
tau2 = times2-times1; %time to travel through cell 2

%sort tau1 and tau2 from low to high and keep their original indices
[tau1_sort,idx1_orig] = sort(tau1);
[tau2_sort,idx2_orig] = sort(tau2);

binsize = length(tau1)/numbins; %approximate number of particles for each bin

%indices and times at which we are separating the particles into bins
idx_cutoff = binsize:binsize:N;
cutoff_time1 = tau1_sort(idx_cutoff); %divide tau1 into equiprobable bins
cutoff_time2 = [0,cutoff_time1]; %use same cutoff times for tau2 as for tau1
cutoff_time2(end) = cutoff_time2(end)+max(tau2); %make sure that last cutoff time is large enough to get max(tau2)
%round cutoff_time2 and tau2_sort to same decimal place for bin search
cutoff_time2 = round(cutoff_time2,4);
tau2_sort = round(tau2_sort,4);

%add index 1 for left edge of first bin
idx_partcutoff1 = [1,idx_cutoff];

binID1_sorted = zeros(N,1); %identify which bin the particles are in (in sorted form tau1_sort)
binID2_sorted = zeros(N,1); %identify which bin the particles are in (in sorted form tau2_sort)

binnumber = 1; %start putting particles in first bin in loop

for jj = 1:numbins
    
    %sort tau1 into equiprobable bins
    binID1_sorted(idx_partcutoff1(jj):idx_partcutoff1(jj+1)-1) = binnumber;
    
    %sort tau2 into the same bins as tau1 
    leftidx_partcutoff2 = min(find(tau2_sort>=cutoff_time2(jj)));
    rightidx_partcutoff2 = max(find(tau2_sort<cutoff_time2(jj+1)));
    binID2_sorted(leftidx_partcutoff2:rightidx_partcutoff2) = binnumber;
    
    if jj == numbins %account for bin number of the last particle
        binID1_sorted(end) = binnumber;
        binID2_sorted(end) = binnumber;
    end
    
    binnumber = binnumber+1;
    
end

%reorganize particles so they are back in their original order
bins1and2 = zeros(N,2); %first column is the particle's bin in cell 1, second column is it's bin in cell 2

for kk = 1:N
    
    %list of bin numbers for each particle in their original order
    bins1and2(idx1_orig(kk),1) = binID1_sorted(kk); %bin for tau1
    bins1and2(idx2_orig(kk),2) = binID2_sorted(kk); %bin for tau2
    
end

bins_tau1 = bins1and2(:,1);
bins_tau2 = bins1and2(:,2);

%%%%% Temporary fix if no tau2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_notau2 = find(bins_tau2==0);
%make bin the same as first bin if no tau 2
bins_tau2(idx_notau2) = bins_tau1(idx_notau2); 
bins1and2(idx_notau2,2) = bins1and2(idx_notau2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop over all the particles and determine which block of the transition
%matrix each tau1/tau2 pair falls into
countTM = zeros(numbins,numbins);
for pp = 1:N 
    
    countTM(bins1and2(pp,1),bins1and2(pp,2))=countTM(bins1and2(pp,1),bins1and2(pp,2))+1;
    
end

getfactor = diff(idx_partcutoff1);
getfactor(1) = getfactor(1)+1; %add one because idx_partcutoff1(1) = 1 not 0

TM = zeros(numbins,numbins);

for mm = 1:numbins
    
    TM(mm,:) = countTM(mm,:)/getfactor(mm);
    
end

%TMnew = flipud(TM');
TMplot = [TM',zeros(numbins,1);zeros(1,numbins+1)];

x = linspace(0,numbins,numbins+1);
y = linspace(0,numbins,numbins+1);
[X,Y] = meshgrid(x,y);

figure(1)
pcolor(X,Y,TMplot)
colorbar
xlabel('Initial state')
ylabel('New state')
%caxis([0 0.15])


idx_state1 = find(bins_tau1==1);
idx_state2 = find(bins_tau1==2);
idx_state3 = find(bins_tau1==3);
idx_state4 = find(bins_tau1==4);
idx_state5 = find(bins_tau1==5);
idx_state6 = find(bins_tau1==6);
idx_state7 = find(bins_tau1==7);
idx_state8 = find(bins_tau1==8);
idx_state9 = find(bins_tau1==9);
idx_state10 = find(bins_tau1==10);
idx_state11 = find(bins_tau1==11);
idx_state12 = find(bins_tau1==12);
idx_state13 = find(bins_tau1==13);
idx_state14 = find(bins_tau1==14);
idx_state15 = find(bins_tau1==15);
idx_state16 = find(bins_tau1==16);
idx_state17 = find(bins_tau1==17);
idx_state18 = find(bins_tau1==18);
idx_state19 = find(bins_tau1==19);
idx_state20 = find(bins_tau1==20);

times_1 = tau1(idx_state1);
times_2 = tau1(idx_state2);
times_3 = tau1(idx_state3);
times_4 = tau1(idx_state4);
times_5 = tau1(idx_state5);
times_6 = tau1(idx_state6);
times_7 = tau1(idx_state7);
times_8 = tau1(idx_state8);
times_9 = tau1(idx_state9);
times_10 = tau1(idx_state10);
times_11 = tau1(idx_state11);
times_12 = tau1(idx_state12);
times_13 = tau1(idx_state13);
times_14 = tau1(idx_state14);
times_15 = tau1(idx_state15);
times_16 = tau1(idx_state16);
times_17 = tau1(idx_state17);
times_18 = tau1(idx_state18);
times_19 = tau1(idx_state19);
times_20 = tau1(idx_state20);

%record initial (y0) and final (y1) y positions across first cell
y0state1 = y0_inlet(idx_state1);
y0state2 = y0_inlet(idx_state2);
y0state3 = y0_inlet(idx_state3);
y0state4 = y0_inlet(idx_state4);
y0state5 = y0_inlet(idx_state5);
y0state6 = y0_inlet(idx_state6);
y0state7 = y0_inlet(idx_state7);
y0state8 = y0_inlet(idx_state8);
y0state9 = y0_inlet(idx_state9);
y0state10 = y0_inlet(idx_state10);
y0state11 = y0_inlet(idx_state11);
y0state12 = y0_inlet(idx_state12);
y0state13 = y0_inlet(idx_state13);
y0state14 = y0_inlet(idx_state14);
y0state15 = y0_inlet(idx_state15);
y0state16 = y0_inlet(idx_state16);
y0state17 = y0_inlet(idx_state17);
y0state18 = y0_inlet(idx_state18);
y0state19 = y0_inlet(idx_state19);
y0state20 = y0_inlet(idx_state20);

[yvec_sort,idxsorty0] = sort(yvec'); %possible y0 values - sort in ascending order for hist
absv_sort = absv(idxsorty0); %corresponding velocity values to y0 - sort in same way as yvec

sum_inv_absv = sum(1./absv_sort); %sum of 1/v1 + 1/v2 + 1/v3 + ... for inverse flux weighting
weightsy0 = 1./absv_sort/sum_inv_absv; %vector = [(1/v1)/sum_inv_absv, (1/v2)/sum_inv_absv, ...]

[y0counts1,y0centers1] = hist(y0state1,yvec_sort);
sum_hist_weights1 = sum(y0counts1.*weightsy0);
Py0_state1 = y0counts1.*weightsy0/sum_hist_weights1; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state1 = cumsum(Py0_state1);

[y0counts2,y0centers2] = hist(y0state2,yvec_sort);
sum_hist_weights2 = sum(y0counts2.*weightsy0);
Py0_state2 = y0counts2.*weightsy0/sum_hist_weights2; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state2 = cumsum(Py0_state2);

[y0counts3,y0centers3] = hist(y0state3,yvec_sort);
sum_hist_weights3 = sum(y0counts3.*weightsy0);
Py0_state3 = y0counts3.*weightsy0/sum_hist_weights3; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state3 = cumsum(Py0_state3);

[y0counts4,y0centers4] = hist(y0state4,yvec_sort);
sum_hist_weights4 = sum(y0counts4.*weightsy0);
Py0_state4 = y0counts4.*weightsy0/sum_hist_weights4; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state4 = cumsum(Py0_state4);

[y0counts5,y0centers5] = hist(y0state5,yvec_sort);
sum_hist_weights5 = sum(y0counts5.*weightsy0);
Py0_state5 = y0counts5.*weightsy0/sum_hist_weights5; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state5 = cumsum(Py0_state5);

[y0counts6,y0centers6] = hist(y0state6,yvec_sort);
sum_hist_weights6 = sum(y0counts6.*weightsy0);
Py0_state6 = y0counts6.*weightsy0/sum_hist_weights6; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state6 = cumsum(Py0_state6);

[y0counts7,y0centers7] = hist(y0state7,yvec_sort);
sum_hist_weights7 = sum(y0counts7.*weightsy0);
Py0_state7 = y0counts7.*weightsy0/sum_hist_weights7; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state7 = cumsum(Py0_state7);

[y0counts8,y0centers8] = hist(y0state8,yvec_sort);
sum_hist_weights8 = sum(y0counts8.*weightsy0);
Py0_state8 = y0counts8.*weightsy0/sum_hist_weights8; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state8 = cumsum(Py0_state8);

[y0counts9,y0centers9] = hist(y0state9,yvec_sort);
sum_hist_weights9 = sum(y0counts9.*weightsy0);
Py0_state9 = y0counts9.*weightsy0/sum_hist_weights9; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state9 = cumsum(Py0_state9);

[y0counts10,y0centers10] = hist(y0state10,yvec_sort);
sum_hist_weights10 = sum(y0counts10.*weightsy0);
Py0_state10 = y0counts10.*weightsy0/sum_hist_weights10; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state10 = cumsum(Py0_state10);

[y0counts11,y0centers11] = hist(y0state11,yvec_sort);
sum_hist_weights11 = sum(y0counts11.*weightsy0);
Py0_state11 = y0counts11.*weightsy0/sum_hist_weights11; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state11 = cumsum(Py0_state11);

[y0counts12,y0centers12] = hist(y0state12,yvec_sort);
sum_hist_weights12 = sum(y0counts12.*weightsy0);
Py0_state12 = y0counts12.*weightsy0/sum_hist_weights12; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state12 = cumsum(Py0_state12);

[y0counts13,y0centers13] = hist(y0state13,yvec_sort);
sum_hist_weights13 = sum(y0counts13.*weightsy0);
Py0_state13 = y0counts13.*weightsy0/sum_hist_weights13; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state13 = cumsum(Py0_state13);

[y0counts14,y0centers14] = hist(y0state14,yvec_sort);
sum_hist_weights14 = sum(y0counts14.*weightsy0);
Py0_state14 = y0counts14.*weightsy0/sum_hist_weights14; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state14 = cumsum(Py0_state14);

[y0counts15,y0centers15] = hist(y0state15,yvec_sort);
sum_hist_weights15 = sum(y0counts15.*weightsy0);
Py0_state15 = y0counts15.*weightsy0/sum_hist_weights15; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state15 = cumsum(Py0_state15);

[y0counts16,y0centers16] = hist(y0state16,yvec_sort);
sum_hist_weights16 = sum(y0counts16.*weightsy0);
Py0_state16 = y0counts16.*weightsy0/sum_hist_weights16; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state16 = cumsum(Py0_state16);

[y0counts17,y0centers17] = hist(y0state17,yvec_sort);
sum_hist_weights17 = sum(y0counts17.*weightsy0);
Py0_state17 = y0counts17.*weightsy0/sum_hist_weights17; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state17 = cumsum(Py0_state17);

[y0counts18,y0centers18] = hist(y0state18,yvec_sort);
sum_hist_weights18 = sum(y0counts18.*weightsy0);
Py0_state18 = y0counts18.*weightsy0/sum_hist_weights18; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state18 = cumsum(Py0_state18);

[y0counts19,y0centers19] = hist(y0state19,yvec_sort);
sum_hist_weights19 = sum(y0counts19.*weightsy0);
Py0_state19 = y0counts19.*weightsy0/sum_hist_weights19; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state19 = cumsum(Py0_state19);

[y0counts20,y0centers20] = hist(y0state20,yvec_sort);
sum_hist_weights20 = sum(y0counts20.*weightsy0);
Py0_state20 = y0counts20.*weightsy0/sum_hist_weights20; %new probability of particle having y_ds=y0 given bin 
cumsumPy0_state20 = cumsum(Py0_state20);

y1state1 = y1_inlet(idx_state1);
y1state2 = y1_inlet(idx_state2);
y1state3 = y1_inlet(idx_state3);
y1state4 = y1_inlet(idx_state4);
y1state5 = y1_inlet(idx_state5);
y1state6 = y1_inlet(idx_state6);
y1state7 = y1_inlet(idx_state7);
y1state8 = y1_inlet(idx_state8);
y1state9 = y1_inlet(idx_state9);
y1state10 = y1_inlet(idx_state10);
y1state11 = y1_inlet(idx_state11);
y1state12 = y1_inlet(idx_state12);
y1state13 = y1_inlet(idx_state13);
y1state14 = y1_inlet(idx_state14);
y1state15 = y1_inlet(idx_state15);
y1state16 = y1_inlet(idx_state16);
y1state17 = y1_inlet(idx_state17);
y1state18 = y1_inlet(idx_state18);
y1state19 = y1_inlet(idx_state19);
y1state20 = y1_inlet(idx_state20);

ntimes1 = length(times_1);
ntimes2 = length(times_2);
ntimes3 = length(times_3);
ntimes4 = length(times_4);
ntimes5 = length(times_5);
ntimes6 = length(times_6);
ntimes7 = length(times_7);
ntimes8 = length(times_8);
ntimes9 = length(times_9);
ntimes10 = length(times_10);
ntimes11 = length(times_11);
ntimes12 = length(times_12);
ntimes13 = length(times_13);
ntimes14 = length(times_14);
ntimes15 = length(times_15);
ntimes16 = length(times_16);
ntimes17 = length(times_17);
ntimes18 = length(times_18);
ntimes19 = length(times_19);
ntimes20 = length(times_20);

save('TMdata.mat')