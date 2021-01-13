function [time,x,timesfunchange,kevalues,indicesfunchange]=robustswingupsol(tspan,x0)
global g l1 l2 m1 m2  flag kp kd ke current previous Ecap Vini T Vdata success in check
options=odeset('Events',@changelyapunov,'AbsTol',1e-13,'RelTol',1e-13);
options2=odeset('Events',@switchtolqr,'AbsTol',1e-13,'RelTol',1e-13);
[time,x,te,xe,ie]=ode15s(@swpencontrol,tspan,x0,options);
timesfunchange=te;
indicesfunchange=length(time);
x0=xe;
tspan=tspan-[0 te];
kevalues=ke;
         t=x0(1);p=x0(2);td=x0(3);pd=x0(4);
      % Setting variables to the right values for the next Lyapunov function
         Ecap=(2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + (2*l1^2*m1*td^2)/3 - g*l1*m1 - 2*g*l1*m2 - g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) - g*l1*m1*cos(t);
         ke=0.98*(3/(4*l2^2*m2))/(abs(Ecap));
        
         Vini=(kp*p^2)/2 + (kd*pd^2)/2 + (ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t))^2)/2;
         flag=0; current=0;previous=1;
         clear -global T Vdata;
in=2;
beginlqrtime=[];
while norm(tspan)>0.1 && isempty(beginlqrtime)  

[timecheck,xcheck,beginlqrtime,beginlqrstate]=ode15s(@swpencontrol,tspan,x0,options2);
% beginlqrtime is the time at which switching to lqr is done considering
% the time of previous function change as time=0
flag=0; current=0;previous=1;

   if isempty(beginlqrtime)

     [time1,x1,te,xe,ie]=ode15s(@swpencontrol,tspan,x0,options);
          if ~(isempty(te))
            timesfunchange(in)=timesfunchange(in-1)+te;
            indicesfunchange(in)=indicesfunchange(in-1)+length(time1);
            x0=xe;
            tspan=tspan-[0 te];
          else
            tspan=tspan-[0 time1(length(time1))];    
          end
     time=cat(1,time,time1+timesfunchange(in-1));
     x=cat(1,x,x1);
     kevalues(in)=ke;
            t=x0(1);p=x0(2);td=x0(3);pd=x0(4);
       % Setting variables to the right values for the next Lyapunov function
            Ecap=(2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + (2*l1^2*m1*td^2)/3 - g*l1*m1 - 2*g*l1*m2 - g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) - g*l1*m1*cos(t);
            ke=0.98*(3/(4*l2^2*m2))/(abs(Ecap));
            Vini=(kp*p^2)/2 + (kd*pd^2)/2 + (ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t))^2)/2;
            flag=0; current=0;previous=1;
            clear -global T Vdata;
     in=in+1;
   else    
     time=cat(1,time,timecheck+timesfunchange(in-1));
     x=cat(1,x,xcheck);
     kevalues(in)=ke;
     x0=beginlqrstate;
     tspan=tspan-[0 beginlqrtime];
     success=1;
     % Adding the last function change time index and value to the respective arrays
     timesfunchange(length(timesfunchange)+1)=time(length(time));
     indicesfunchange(length(indicesfunchange)+1)=length(time);
   end   
end

if success==1
    optionslqr=odeset('AbsTol',1e-13,'RelTol',1e-13);
    [timelqr,xlqr]=ode15s('swpencontrol',tspan,x0,optionslqr);
    time=cat(1,time,timelqr+time(length(time)));
    x=cat(1,x,xlqr);
end
