%reflects particles off top and bottom boundary
    
    %top boundary
    out=find(yA>Ly);    
    yA(out)=2*Ly-yA(out); 
    clear out 
    
    out=find(yB>Ly);    
    yB(out)=2*Ly-yB(out);
    clear out 
    

    %bottom boundary
    out=find(yA<1);    
    yA(out)=1-yA(out);
    clear out

    out=find(yB<1);    
    yB(out)=2-yB(out);
    clear out
