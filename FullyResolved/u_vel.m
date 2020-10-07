%rk4.m

function [u] = u_vel(x,y)

global V

xx=floor(x);
yy=floor(y);

%xxx=ceil(x);
xxx=xx + 1;
%yyy=ceil(y);
yyy=yy + 1;

u=zeros(size(xx));

for ii=1:length(xx)
    
        u(ii)=V.x(xx(ii),yy(ii))+(x(ii)-xx(ii))./(xxx(ii)-xx(ii)).*(V.x(xxx(ii),yy(ii))-V.x(xx(ii),yy(ii)));
    
end
	