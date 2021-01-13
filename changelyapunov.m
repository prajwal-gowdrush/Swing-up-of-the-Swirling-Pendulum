function [value,isterminal,direction] = changelyapunov(time,x)
global Vdata T current previous Vini flag ke Ecap l2 m2 in check
t=x(1);
p=x(2);
td=x(3);
pd=x(4);
p=2*pi*(p/(2*pi)-fix(p/(2*pi)));p=pi*(p/pi-2*fix(p/pi));% setting the range of p to [-pi,pi]
t=2*pi*(t/(2*pi)-floor(t/(2*pi)));
t_err=t-pi;
if flag==1
value =((Vdata(previous)-Vdata(current))/(T(current)-T(previous)))-(0.1*(0.9*Vini)/(T(1)));
else
value=1; 
end
isterminal = 1;% Stop the integration
direction = 0;   % Any Direction