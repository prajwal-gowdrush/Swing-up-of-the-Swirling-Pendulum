function xdot=swpencontrol(time,x)
global g l1 l2 m1 m2 u C flag kp kd ke effort current previous Vini T Vdata Ecap success
t=x(1);
p=x(2);
td=x(3);
pd=x(4);
p=2*pi*(p/(2*pi)-fix(p/(2*pi)));p=pi*(p/pi-2*fix(p/pi));% setting the range of p to [-pi,pi]
t=2*pi*(t/(2*pi)-floor(t/(2*pi)));
t_err=t-pi;
distnorm=sqrt((t_err)^2+p^2+td^2+pd^2);

V=(kp*p^2)/2 + (kd*pd^2)/2 + (ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t))^2)/2;
if V<(0.1*Vini)
current=current+1;
T(current)=time;
Vdata(current)=V;
     while (T(current)-T(previous))>2
          previous=previous+1;
          flag=1;
     end        
end

% Ecap=(2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + (2*l1^2*m1*td^2)/3 - g*l1*m1 - 2*g*l1*m2 - g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) - g*l1*m1*cos(t);
%for finding optimal ke values

xdot(1,1)=x(3);
xdot(2,1)=x(4);
%For uncontrolled system
% xdot(3,1)=-(16*m2*sin(p)*l2^2*pd*td*cos(p) - 12*l1*m2*sin(p)*l2*pd^2 + 12*l1*m2*sin(p)*l2*td^2*cos(p)^2 + 6*g*m2*cos(t)*sin(p)*l2 - 9*g*l1*m2*sin(t)*cos(p)^2 + 6*g*l1*m1*sin(t) + 12*g*l1*m2*sin(t))/(8*l1^2*m1 + 24*l1^2*m2 + 8*l2^2*m2 - 18*l1^2*m2*cos(p)^2 - 8*l2^2*m2*cos(p)^2);
% xdot(4,1)=(6*g*l2^2*m2*cos(p)^3*sin(t) + 8*l2^3*m2*td^2*cos(p)*sin(p) - 8*l2^3*m2*td^2*cos(p)^3*sin(p) + 3*g*l1^2*m1*cos(p)*sin(t) - 6*g*l2^2*m2*cos(p)*sin(t) - 18*l1^2*l2*m2*pd^2*cos(p)*sin(p) + 8*l1^2*l2*m1*td^2*cos(p)*sin(p) + 24*l1^2*l2*m2*td^2*cos(p)*sin(p) + 9*g*l1*l2*m2*cos(p)*cos(t)*sin(p) + 24*l1*l2^2*m2*pd*td*cos(p)^2*sin(p))/(8*l2^3*m2 + 8*l1^2*l2*m1 + 24*l1^2*l2*m2 - 8*l2^3*m2*cos(p)^2 - 18*l1^2*l2*m2*cos(p)^2);

%For controlled system
% if flag==0 &&
% sqrt((t-pi)^2+(C(2)/C(1)*p)^2+(C(3)/C(1)*td)^2+(C(4)/C(1)*pd)^2)>=ellipdisteq
% %Closed Loop equations for swing-up control
% xdot(3,1)=-((9*l1*cos(p)*(pd + kp*p + (kd*((2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p)))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2))))/(ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(3*l1^2*td^2 + 3*l1*l2*pd*td*cos(p) + l2^2*pd^2 - l2^2*td^2*cos(p)^2 + l2^2*td^2))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t)) - kd/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2))) - 12*l1*l2^2*m2*pd^2*sin(p) + 6*g*l2^2*m2*cos(t)*sin(p) + 6*g*l1*l2*m1*sin(t) + 12*g*l1*l2*m2*sin(t) - 9*g*l1*l2*m2*cos(p)^2*sin(t) + 16*l2^3*m2*pd*td*cos(p)*sin(p) + 12*l1*l2^2*m2*td^2*cos(p)^2*sin(p))/(8*l2^3*m2 + 8*l1^2*l2*m1 + 24*l1^2*l2*m2 - 8*l2^3*m2*cos(p)^2 - 18*l1^2*l2*m2*cos(p)^2);
% xdot(4,1)=((pd + kp*p + (kd*((2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p)))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2)))/(ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(3*l1^2*td^2 + 3*l1*l2*pd*td*cos(p) + l2^2*pd^2 - l2^2*td^2*cos(p)^2 + l2^2*td^2))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t)) - kd/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2))) + (2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2));
% else
% flag=1;    
% u=-C*(x-[pi;0;0;0]); % For LQR stabilisation about upright equilibrium
% xdot(3,1)=-(9*l1*u*cos(p) - 12*l1*l2^2*m2*pd^2*sin(p) + 6*g*l2^2*m2*cos(t)*sin(p) + 6*g*l1*l2*m1*sin(t) + 12*g*l1*l2*m2*sin(t) - 9*g*l1*l2*m2*cos(p)^2*sin(t) + 16*l2^3*m2*pd*td*cos(p)*sin(p) + 12*l1*l2^2*m2*td^2*cos(p)^2*sin(p))/(8*l2^3*m2 + 8*l1^2*l2*m1 + 24*l1^2*l2*m2 - 8*l2^3*m2*cos(p)^2 - 18*l1^2*l2*m2*cos(p)^2);
% xdot(4,1)=(u + (2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2));
% end
%% Closed Loop equations for LQR (with c1 c2 c3 c4 specified)
% c1=C(1); c2=C(2); c3=C(3); c4=C(4);
% xdot(3,1)=-(6*g*l2^2*m2*cos(t)*sin(p) - 12*l1*l2^2*m2*pd^2*sin(p) - 9*l1*cos(p)*(c2*p + c4*pd + c3*td - c1*(pi - t)) + 6*g*l1*l2*m1*sin(t) + 12*g*l1*l2*m2*sin(t) - 9*g*l1*l2*m2*cos(p)^2*sin(t) + 16*l2^3*m2*pd*td*cos(p)*sin(p) + 12*l1*l2^2*m2*td^2*cos(p)^2*sin(p))/(8*l2^3*m2 + 8*l1^2*l2*m1 + 24*l1^2*l2*m2 - 8*l2^3*m2*cos(p)^2 - 18*l1^2*l2*m2*cos(p)^2);
% xdot(4,1)=-(c2*p + c4*pd + c3*td - c1*(pi - t) - (2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 + g*l2*m2*cos(p)*sin(t) - (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) - 2*l1*l2*m2*pd*td*sin(p))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2));
%% For controlled system 
if success==1
     u=-C*([t;p;td;pd]-[pi;0;0;0]);  %For LQR stabilisation about upright equilibrium
         if u>1000 
               disp('Disturbance is too large to stabilise');
           return 
         end
xdot(3,1)=-(9*l1*u*cos(p) - 12*l1*l2^2*m2*pd^2*sin(p) + 6*g*l2^2*m2*cos(t)*sin(p) + 6*g*l1*l2*m1*sin(t) + 12*g*l1*l2*m2*sin(t) - 9*g*l1*l2*m2*cos(p)^2*sin(t) + 16*l2^3*m2*pd*td*cos(p)*sin(p) + 12*l1*l2^2*m2*td^2*cos(p)^2*sin(p))/(8*l2^3*m2 + 8*l1^2*l2*m1 + 24*l1^2*l2*m2 - 8*l2^3*m2*cos(p)^2 - 18*l1^2*l2*m2*cos(p)^2);
xdot(4,1)=(u + (2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2));
else
u=(pd + kp*p + (kd*((2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p)))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2)))/(ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(3*l1^2*td^2 + 3*l1*l2*pd*td*cos(p) + l2^2*pd^2 - l2^2*td^2*cos(p)^2 + l2^2*td^2))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t)) - kd/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2)));
%Above are the control inputs for swing-up   
xdot(3,1)=-(9*l1*u*cos(p) - 12*l1*l2^2*m2*pd^2*sin(p) + 6*g*l2^2*m2*cos(t)*sin(p) + 6*g*l1*l2*m1*sin(t) + 12*g*l1*l2*m2*sin(t) - 9*g*l1*l2*m2*cos(p)^2*sin(t) + 16*l2^3*m2*pd*td*cos(p)*sin(p) + 12*l1*l2^2*m2*td^2*cos(p)^2*sin(p))/(8*l2^3*m2 + 8*l1^2*l2*m1 + 24*l1^2*l2*m2 - 8*l2^3*m2*cos(p)^2 - 18*l1^2*l2*m2*cos(p)^2);
xdot(4,1)=(u + (2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2));

%Closed-Loop equations for swing-up control    
% xdot(3,1)=-((9*l1*cos(p)*(pd + kp*p + (kd*((2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p)))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2))))/(ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(3*l1^2*td^2 + 3*l1*l2*pd*td*cos(p) + l2^2*pd^2 - l2^2*td^2*cos(p)^2 + l2^2*td^2))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t)) - kd/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2))) - 12*l1*l2^2*m2*pd^2*sin(p) + 6*g*l2^2*m2*cos(t)*sin(p) + 6*g*l1*l2*m1*sin(t) + 12*g*l1*l2*m2*sin(t) - 9*g*l1*l2*m2*cos(p)^2*sin(t) + 16*l2^3*m2*pd*td*cos(p)*sin(p) + 12*l1*l2^2*m2*td^2*cos(p)^2*sin(p))/(8*l2^3*m2 + 8*l1^2*l2*m1 + 24*l1^2*l2*m2 - 8*l2^3*m2*cos(p)^2 - 18*l1^2*l2*m2*cos(p)^2);
% xdot(4,1)=((pd + kp*p + (kd*((2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p)))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2)))/(ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(3*l1^2*td^2 + 3*l1*l2*pd*td*cos(p) + l2^2*pd^2 - l2^2*td^2*cos(p)^2 + l2^2*td^2))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t)) - kd/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2))) + (2*m2*(2*l2^2*td^2*cos(p)*sin(p) - 3*l1*l2*pd*td*sin(p)))/3 - g*l2*m2*cos(p)*sin(t) + (2*l1*l2*m2*cos(p)*(8*m2*td*cos(p)*sin(p)*l2^2*pd - 6*l1*m2*sin(p)*l2*pd^2 + 3*g*m2*cos(t)*sin(p)*l2 + 3*g*l1*m1*sin(t) + 6*g*l1*m2*sin(t)))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2) + 2*l1*l2*m2*pd*td*sin(p))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(p)^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(p)^2));
end

