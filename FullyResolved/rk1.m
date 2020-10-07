%rk4.m

function [xnew, ynew] = rk1 (x0, y0, dt);
global V D

    xe1=sqrt(2*D*dt)*randn(size(x0));
    ye1=sqrt(2*D*dt)*randn(size(y0));

    x1 = x0+dt*u_vel(x0,y0)+xe1;
    y1 = y0+dt*v_vel(x0,y0)+ye1;
  
    xnew = x1;%x0+dt*a31*u_vel(x0,y0)+dt*a32*u_vel(x1,y1)+dt*a33*u_vel(x2,y2)+b31*xe1+b32*xe2+b33*xe3;
    ynew = y1;%y0+dt*a31*v_vel(x0,y0)+dt*a32*v_vel(x1,y1)+dt*a33*v_vel(x2,y2)+b31*ye1+b32*ye2+b33*ye3;
