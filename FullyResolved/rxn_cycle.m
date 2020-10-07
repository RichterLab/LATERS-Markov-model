%reaction step

xdead = zeros(size(xA));   %tracker to kill of A particles that have reacted
reactedpart = zeros(length(IdxTot),3);   %1st column = x position of reacted particles
                                      %2nd column = y position of reacted particles
                                      %3rd column = kk (time step in which reaction occurred)
react_velA = zeros(length(IdxTot),2);
react_velB = zeros(length(IdxTot),2);
recentloc_xA = zeros(length(IdxTot),10);
recentloc_yA = zeros(length(IdxTot),10);
recentloc_xB = zeros(length(IdxTot),10);
recentloc_yB = zeros(length(IdxTot),10);
    
for jj=1:length(IdxTot) %loop over A particles to determine which (if any) B particles they react with
   
    Bidx = IdxTot{jj}; %vector of the indices of B particles within range of the jj-th A particle
    
        xBtemp = xB(Bidx); %locations of B points within range of the particular A point
        yBtemp = yB(Bidx);
       
        X = xA(jj)*ones(1,length(Bidx)); %the position of the A particle we are looking at
        Y = yA(jj)*ones(1,length(Bidx));
        
        %calculate minimum distance between this A and all B particles accounting for periodic
        %domain.
        s = sqrt((X-xBtemp).^2+(Y-yBtemp).^2);
        s = min(s);
        
        %calculate the probability of reaction of each particle pair given
        %the distance between them
        Pf = Pr*1/(8*pi*D*dt)*exp(-s.^2/(8*D*dt));

        %generate a random number the size of each particle pair
        RP = Pf-rand(size(Pf));

        %identify the most probable of all possible reactions
        idxreact_temp = find(RP==max(RP)); %index from list of B particles within range of A particle 

        %kill the B particle that partook in the most probable reaction
        if RP(idxreact_temp)>0
            
            xdead(jj) = 1;
            idxreact = Bidx(idxreact_temp); %index of B particle that reacts
            
            %record location of reaction (occurs at midpoint between the
            %two reacting particles), velocity of reacting particles, and
            %the 10 most recent locations of the particles
            reactedpart(jj,1) = (xA(jj)+xB(idxreact))/2; %x position of reaction
            reactedpart(jj,2) = (yA(jj)+yB(idxreact))/2; %y position of reaction
            reactedpart(jj,3) = kk;
            
            recentloc_xA(jj,:) = xArecent(jj,:);
            recentloc_yA(jj,:) = yArecent(jj,:);
            recentloc_xB(jj,:) = xBrecent(idxreact,:);
            recentloc_yB(jj,:) = yBrecent(idxreact,:);
            
            %mark B particle for reaction
            xB(idxreact) = NaN;
            yB(idxreact) = NaN;
            
            %create new C particle
            xC = [xC,reactedpart(jj,1)];
            yC = [yC,reactedpart(jj,2)];
            
        end
        
        %If they don't react, nothing happens. The particles stay in the
        %system and move by random walk in the next time step.
    
end
    
    %kill all particles that took place in a reaction
    notdead = find(xdead==0); %indices of A particles that haven't reacted
    dead = find(xdead==1);    %indices of A particles that have reacted
    
    %update the total count of particles that have reacted with the number 
    %that have reacted during this timestep
    countprevious = countreact;
    countreact = countreact+length(dead);
    
    %store locations of reactions, velocities, and recent locations of
    %reacting particles
    react_loc((countprevious+1):countreact,:) = reactedpart(dead,:);
    xArecent_loc((countprevious+1):countreact,:) = recentloc_xA(dead,:);
    yArecent_loc((countprevious+1):countreact,:) = recentloc_yA(dead,:);
    xBrecent_loc((countprevious+1):countreact,:) = recentloc_xB(dead,:);
    yBrecent_loc((countprevious+1):countreact,:) = recentloc_yB(dead,:);
    
    xA = xA(notdead);
    yA = yA(notdead);
    xB = xB(~isnan(xB));
    yB = yB(~isnan(yB));
    
    xArecent = xArecent(notdead,:);
    yArecent = yArecent(notdead,:);
    xBrecent = xBrecent(~isnan(xB),:);
    yBrecent = yBrecent(~isnan(yB),:);
   