%rk4.m

function [v] = v_vel(x,y)

global V

xx=floor(x);
yy=floor(y);

%xxx=ceil(x);
xxx=xx + 1;
%yyy=ceil(y);
yyy=yy + 1;

v=zeros(size(xx));

for ii=1:length(xx)
   
    v(ii)=V.y(xx(ii),yy(ii))+(y(ii)-yy(ii))./(yyy(ii)-yy(ii)).*(V.y(xx(ii),yyy(ii))-V.y(xx(ii),yy(ii)));
   
end
    	