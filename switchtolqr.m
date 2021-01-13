function [value,isterminal,direction] = switchtolqr(time,x)
global ke
t=x(1);
p=x(2);
td=x(3);
pd=x(4);
p=2*pi*(p/(2*pi)-fix(p/(2*pi)));p=pi*(p/pi-2*fix(p/pi));% setting the range of p to [-pi,pi]
t=2*pi*(t/(2*pi)-floor(t/(2*pi)));
t_err=t-pi;
distnorm=sqrt((t_err)^2+p^2+td^2+pd^2);
value=distnorm-0.0039;
isterminal = 1;   % Stop the integration
direction = 0;   % Any Direction