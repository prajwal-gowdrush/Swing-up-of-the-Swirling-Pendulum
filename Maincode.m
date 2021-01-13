 %Main Code - Upright equilibrium
 clear all;
 close all;
 clc;
global g l1 l2 m1 m2 u C flag kp kd ke effort current previous Ecap Vini T Vdata success check
g=9.81;
flag=0;
current=0;
previous=1;
success=0;
check=0;
%% Input the parameter values
% choice = menu('Assign default link lengths and masses?','Yes','No');
choice=1;
if choice==1
    l1=0.2;
    l2=0.4;
    m1=0.2;
    m2=0.4;
elseif choice==2
    l1=input('The length of link 1(in metres)= ')/2;
    l2=input('The length of link 2(in metres)= ')/2;
    m1=input('The mass of link 1(in kg)=');
    m2=input('The mass of link 2(in kg)='); 
end
%% Input the initial conditions
% x0(1,1)=input('Initial angular position of link 1= ');
% x0(2,1)=input('Initial angular position of link 2= ');
% x0(3,1)=input('Initial angular velocity of link 1= ');
% x0(4,1)=input('Initial angular velocity of link 2= ');
% simtim=input('Simulation is needed for a time of ___seconds ');
% tspan=[0 simtim];

% x0(1,1)=0.1395+0.01;
% x0(2,1)=-2.9652-0.01;

x0(1,1)=pi/4;
x0(2,1)=0;

x0(3,1)=0;
x0(4,1)=0;
simtim=50;
tspan=[0 simtim];

%% Initialisation of some variables for swing-up control
t=x0(1);p=x0(2);td=x0(3);pd=x0(4);
Ecap=(2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + (2*l1^2*m1*td^2)/3 - g*l1*m1 - 2*g*l1*m2 - g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) - g*l1*m1*cos(t);
kp=1;kd=1;
% t=acos((-981/500+(0.95*(3/(4*l2^2*m2)))/ke)*500/981);%calculating theta from ke
ke=0.95*(3/(4*l2^2*m2))/(abs(Ecap));
% ke=0.95;
Vini=(kp*p^2)/2 + (kd*pd^2)/2 + (ke*(g*l1*m1 - (2*l1^2*m1*td^2)/3 - (2*m2*(l2^2*pd^2 + 3*l1^2*td^2 + l2^2*td^2 - l2^2*td^2*cos(p)^2 + 3*l1*l2*pd*td*cos(p)))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(t) - l2*sin(p)*sin(t)) + g*l1*m1*cos(t))^2)/2;
%% Calculations
%LQR Matrices
% u=0;
A=[0 0 1 0;0 0 0 1;(l1*(6*g*m1 + 3*g*m2))/(8*l1^2*m1 + 6*l1^2*m2) (3*g*l2*m2)/(l1^2*(4*m1 + 3*m2)) 0 0;-(3*g*m1)/(2*l2*(4*m1 + 3*m2)) -(9*g*m2)/(8*l1*m1 + 6*l1*m2) 0 0];
B=[0;0;-(9*l1)/(8*l1^2*l2*m1+6*l1^2*l2*m2); 1/((4*l2^2*m2)/3-(12*l1^2*l2^2*m2^2)/(4*l1^2*m1+12*l1^2*m2))];
C=lqr(A,B,1*eye(4),1,zeros(4,1));
%% Solution to the differential equations

% Robust solution for complete swing-up 
[time,x,timesfunchange,kevalues,indicesfunchange]=robustswingupsol(tspan,x0);

% % A single step solution for a constant ke
% options = odeset('AbsTol',1e-13,'RelTol',1e-13);
% [time,x]=ode15s('swpencontrol',tspan,x0,options);

% Extracting information from the solution matrix of the ode
theta=x(:,1);phi=x(:,2);thetadot=x(:,3);phidot=x(:,4);
readtheta=2*pi*(theta/(2*pi)-floor(theta/(2*pi)));
thetaerr=readtheta-pi;
phierr=2*pi*(phi/(2*pi)-fix(phi/(2*pi)));phierr=pi*(phierr/pi-2*fix(phierr/pi));

%link coordinates calculations for plotting
 x1=2*l1*sin(theta);
 y1=0*theta;
 z1=-2*l1*cos(theta);
 x2=2*l2*sin(phi).*cos(theta)+x1;
 y2=2*l2*cos(phi);
 z2=2*l2*sin(phi).*sin(theta)+z1;

 n=size(time);
 n=n(1);wait=(simtim/n);
 

%  
%  Settling Time calculation
%   compvec=abs(pi-theta)>0*theta+abs(0.05*(pi-x0(1)));
%  for i=1:n
%   if sum(compvec(i:n))==0
%          settlingtime=time(i);
%          break;
%   end
%  end
%  
%% MODEL VALIDATION (energy conservation must hold when no control is appied)
 energy=(2*m2*(3*l1^2.*thetadot.^2 + 3*l1*l2*phidot.*thetadot.*cos(phi) + l2^2.*phidot.^2 - l2^2*thetadot.^2.*cos(phi).^2 + l2^2*thetadot.^2))/3 + (2*l1^2*m1*thetadot.^2)/3 - g*m2*(2*l1*cos(theta) - l2*sin(phi).*sin(theta)) - g*l1*m1*cos(theta);
 kenergy=(2*m2*(l2^2*phidot.^2 + 3*l1^2*thetadot.^2 + l2^2*thetadot.^2 - l2^2*thetadot.^2.*cos(phi).^2 + 3*l1*l2*phidot.*thetadot.*cos(phi)))/3 + (2*l1^2*m1*thetadot.^2)/3;
%  deviation=energy-energy(1);
%  deviation=sum(abs(deviation))/n;
%  avg_percent_error_energy=deviation/energy(1)*100
%  figure;
%  plot(time,energy);
%  title('Energy of the free system');
%  ylabel('Energy');
%  xlabel('Time');
%  axis([0 simtim energy(1)-2 energy(1)+2]);

%% Plot of Joint-Angles vs Time for swing-up + capture
linethick=1.5;%used for all plots

figure;
plot(time,theta,'r',time,phi,'b','LineWidth',linethick);
 xlabel('Time(in seconds)');
 ylabel('Joint angles(in radians)');
%  title('Swing-Up');
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
hold on;
 for verline=1:length(timesfunchange)
  plot([timesfunchange(verline) timesfunchange(verline)],get(ax,'YLim'),'--m','LineWidth',1.4);
 end
 legstr=strcat('ke=',num2str(kevalues(1)),', ke=',num2str(kevalues(2)),', ke=',num2str(kevalues(3)), ', LQR');
 legend({'\theta','\phi',legstr},'FontSize',16)
 
 %% Plot of Link Velocities
 figure;
plot(time,thetadot,'r',time,phidot,'b','LineWidth',linethick);
 xlabel('Time(in seconds)');
 ylabel('Joint Velocities(in radians/second)');
%  title('Swing-Up');
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
hold on;
 for verline=1:length(timesfunchange)
  plot([timesfunchange(verline) timesfunchange(verline)],get(ax,'YLim'),'--m','LineWidth',linethick);
 end
 legstr=strcat('ke=',num2str(kevalues(1)),', ke=',num2str(kevalues(2)),', ke=',num2str(kevalues(3)), ', LQR');
 legend({'\theta_{dot}','\phi_{dot}',legstr},'FontSize',16)
%temporary
% plot(time,theta,'r');
%  xlabel('Time');
%  ylabel('Theta');
% %  title('Check');
%  legend({'Theta'},'FontSize',16)
% grid on;
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 1;
% ax.FontSize = 16;
% ax.LineWidth = 1.4;
% 
% figure;
% plot(time,phi,'b');
%  xlabel('Time');
%  ylabel('Phi');
% %  title('Swing-Up');
%  legend({'Phi'},'FontSize',16)
% grid on;
% ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 1;
% ax.FontSize = 16;
% ax.LineWidth = 1.4;

 %for swing-up control
phi_inrange=2*pi*(phi/(2*pi)-fix(phi/(2*pi)));
phi_inrange=pi*(phi_inrange/pi-2*fix(phi_inrange/pi));

%% for swing-up control (temp) 
% phi_inrange=2*pi*(phi/(2*pi)-fix(phi/(2*pi)));
% phi_inrange=pi*(phi_inrange/pi-2*fix(phi_inrange/pi));
% %Lyapunov function calculation
% V=0.5*(kp*phi_inrange.^2+kd*phidot.^2+ke*(energy-g*l1*m1-g*m2*(2*l1 + l2)).^2);
%  figure;
%  plot(time,V); 
%  title('Lyapunov function');
%  ylabel('V');
%  xlabel('Time');
%  axis([0 simtim 0 1.5*max(V)]);
% 
%% Plots of variation of Lyapunov functions with time
for i=1:(length(indicesfunchange))
   % %Lyapunov function calculation
   if i>1
phi_inrange_int=phi_inrange(indicesfunchange(i-1):indicesfunchange(i));
phidot_int=phidot(indicesfunchange(i-1):indicesfunchange(i));
thetadot_int=thetadot(indicesfunchange(i-1):indicesfunchange(i));
theta_int=theta(indicesfunchange(i-1):indicesfunchange(i));
energy_int=energy(indicesfunchange(i-1):indicesfunchange(i));
time_int=time(indicesfunchange(i-1):indicesfunchange(i));
   else
phi_inrange_int=phi_inrange(1:indicesfunchange(i));
phidot_int=phidot(1:indicesfunchange(i));
thetadot_int=thetadot(1:indicesfunchange(i));
theta_int=theta(1:indicesfunchange(i));
energy_int=energy(1:indicesfunchange(i));
time_int=time(1:indicesfunchange(i));
   end
V=0.5*(kp*phi_inrange_int.^2+kd*phidot_int.^2+kevalues(i)*(energy_int-g*l1*m1-2*g*l1*m2).^2);
 figure;
 plot(time_int,V,'b','LineWidth',1.4); 
%  title('Lyapunov function');
 ylabel('Lyapunov Function(V)');
 xlabel('Time(in seconds)');
 if i>1
 axis([timesfunchange(i-1) timesfunchange(i) 0 1.5*max(V)]);
 else
 axis([0 timesfunchange(i) 0 1.5*max(V)]);
 end
 grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
legstr2=strcat('kp=1, kd=1, ke=', num2str(kevalues(i)));
legend({legstr2},'FontSize',16);
% Effort calculation
 effort1=zeros(length(time_int),1);
 for ie=1:length(time_int)
effort1(ie)=(phidot_int(ie) + kp*phi_inrange_int(ie) + (kd*((2*m2*(2*l2^2*thetadot_int(ie)^2*cos(phi_inrange_int(ie))*sin(phi_inrange_int(ie)) - 3*l1*l2*phidot_int(ie)*thetadot_int(ie)*sin(phi_inrange_int(ie))))/3 - g*l2*m2*cos(phi_inrange_int(ie))*sin(theta_int(ie)) + (2*l1*l2*m2*cos(phi_inrange_int(ie))*(8*m2*thetadot_int(ie)*cos(phi_inrange_int(ie))*sin(phi_inrange_int(ie))*l2^2*phidot_int(ie) - 6*l1*m2*sin(phi_inrange_int(ie))*l2*phidot_int(ie)^2 + 3*g*m2*cos(theta_int(ie))*sin(phi_inrange_int(ie))*l2 + 3*g*l1*m1*sin(theta_int(ie)) + 6*g*l1*m2*sin(theta_int(ie))))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi_inrange_int(ie))^2) + 2*l1*l2*m2*phidot_int(ie)*thetadot_int(ie)*sin(phi_inrange_int(ie))))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(phi_inrange_int(ie))^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi_inrange_int(ie))^2)))/(kevalues(i)*(g*l1*m1 - (2*l1^2*m1*thetadot_int(ie)^2)/3 - (2*m2*(3*l1^2*thetadot_int(ie)^2 + 3*l1*l2*phidot_int(ie)*thetadot_int(ie)*cos(phi_inrange_int(ie)) + l2^2*phidot_int(ie)^2 - l2^2*thetadot_int(ie)^2*cos(phi_inrange_int(ie))^2 + l2^2*thetadot_int(ie)^2))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(theta_int(ie)) - l2*sin(phi_inrange_int(ie))*sin(theta_int(ie))) + g*l1*m1*cos(theta_int(ie))) - kd/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(phi_inrange_int(ie))^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi_inrange_int(ie))^2)));
 % effort(ie)=(phidot(ie) + kp*phi(ie) + (kd*((2*m2*(2*l2^2*thetadot(ie)^2*cos(phi(ie))*sin(phi(ie)) - 3*l1*l2*phidot(ie)*thetadot(ie)*sin(phi(ie))))/3 - g*l2*m2*cos(phi(ie))*sin(theta(ie)) + (2*l1*l2*m2*cos(phi(ie))*(8*m2*thetadot(ie)*cos(phi(ie))*sin(phi(ie))*l2^2*phidot(ie) - 6*l1*m2*sin(phi(ie))*l2*phidot(ie)^2 + 3*g*m2*cos(theta(ie))*sin(phi(ie))*l2 + 3*g*l1*m1*sin(theta(ie)) + 6*g*l1*m2*sin(theta(ie))))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi(ie))^2) + 2*l1*l2*m2*phidot(ie)*thetadot(ie)*sin(phi(ie))))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(phi(ie))^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi(ie))^2)))/(ke*(g*l1*m1 - (2*l1^2*m1*thetadot(ie)^2)/3 - (2*m2*(3*l1^2*thetadot(ie)^2 + 3*l1*l2*phidot(ie)*thetadot(ie)*cos(phi(ie)) + l2^2*phidot(ie)^2 - l2^2*thetadot(ie)^2*cos(phi(ie))^2 + l2^2*thetadot(ie)^2))/3 + 2*g*l1*m2 + g*m2*(2*l1*cos(theta(ie)) - l2*sin(phi(ie))*sin(theta(ie))) + g*l1*m1*cos(theta(ie))) - kd/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(phi(ie))^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi(ie))^2)));
 end 
 effort=cat(1,effort,effort1);
end

%  Control effort calculation for LQR
a=indicesfunchange(length(indicesfunchange));
b=length(thetaerr);
x_int=[(thetaerr(a:b))';(phierr(a:b))';(thetadot(a:b))';(phidot(a:b))'];
effort1=(-C*(x_int))';
effort=cat(1,effort,effort1);
effort=effort(1:size(time,1),1);
 %Vdot calculation
  phidoubledot=zeros(size(time,1),1);
 for ie=1:size(time,1)
   phidoubledot(ie)=(effort(ie) + (2*m2*(2*l2^2*thetadot(ie)^2*cos(phi(ie))*sin(phi(ie)) - 3*l1*l2*phidot(ie)*thetadot(ie)*sin(phi(ie))))/3 - g*l2*m2*cos(phi(ie))*sin(theta(ie)) + (2*l1*l2*m2*cos(phi(ie))*(8*m2*thetadot(ie)*cos(phi(ie))*sin(phi(ie))*l2^2*phidot(ie) - 6*l1*m2*sin(phi(ie))*l2*phidot(ie)^2 + 3*g*m2*cos(theta(ie))*sin(phi(ie))*l2 + 3*g*l1*m1*sin(theta(ie)) + 6*g*l1*m2*sin(theta(ie))))/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi(ie))^2) + 2*l1*l2*m2*phidot(ie)*thetadot(ie)*sin(phi(ie)))/((4*l2^2*m2)/3 - (12*l1^2*l2^2*m2^2*cos(phi(ie))^2)/(4*l1^2*m1 + 12*l1^2*m2 + 4*l2^2*m2 - 4*l2^2*m2*cos(phi(ie))^2));
 end 
%  Vdot=phidot.*(kp*phi+kd*phidoubledot+ke*(energy-g*l1*m1-2*g*l1*m2).*effort);
Vdot=-phidot.*phidot;
 %% Plots of Stabilisation about the upward equilibrium
%  subplot(1,2,1)
%  plot(time,theta,'b');
%  axis([0 simtim min([pi pi+1.5*min(theta-pi)]) pi+1.5*max(theta-pi)]);
%  xlabel('Time');
%  ylabel('Theta');
%  title('Link 1');
%   
%  subplot(1,2,2)
%  plot(time,phi,'r');
%  axis([0 simtim min([0 1.5*min(phi)]) 1.5*max(phi)]);
%  xlabel('Time');
%  ylabel('Phi');
%  title('Link 2');

 %% In a single plot
% plot(time,theta-pi,'r',time,phi,'b','LineWidth',1.5);
%  xlabel('Time');
%  ylabel('Joint angles(rad)');
% %  title('LQR Stabilisation');
%  legend({'Disturbance in Theta','Disturbance in Phi'},'FontSize',16);
%  grid on;
%  ax = gca;
% ax.GridLineStyle = ':';
% ax.GridAlpha = 1;
% ax.FontSize = 16;
% ax.LineWidth = 1.4;

%temporary lines
% x_temp=[(thetaerr)';(phierr)';(thetadot)';(phidot)'];
% effort=(-C*(x_temp))';

% % Control effort plot
 figure;
 plot(time,effort(1:length(time)),'b','LineWidth',linethick);
%  axis([0 simtim min([0 min(effort)]) max(effort)]);
 xlabel('Time(in seconds)');
 ylabel('Control Effort(in N-m)');
%  title('Control Effort');
 grid on;
 ax = gca;
ax.GridLineStyle = ':';  
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
hold on;
 for verline=1:length(timesfunchange)
  plot([timesfunchange(verline) timesfunchange(verline)],get(ax,'YLim'),'--m','LineWidth',linethick);
 end
 legend({'u(t)',legstr},'FontSize',16)
% %  
% %Energy plot
figure;
 plot(time,energy,'b','LineWidth',linethick);
%  axis([0 simtim min([0 min(effort)]) max(effort)]);
 xlabel('Time(in seconds)');
 ylabel('Total Energy(in Joules)');
%  title('Energy');
 grid on;
 ax = gca;
ax.GridLineStyle = ':';  
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
hold on;
 for verline=1:length(timesfunchange)
  plot([timesfunchange(verline) timesfunchange(verline)],get(ax,'YLim'),'--m','LineWidth',linethick);
 end
 legend({'E(t)',legstr},'FontSize',16)

 %% Vdot plot
  figure;
 plot(time,Vdot,'b','LineWidth',linethick);
 axis([0 simtim min([0 1.25*min(Vdot)]) 1.25*max(Vdot)]);
 xlabel('Time');
 ylabel('Vdot');
  grid on;
 ax = gca;
ax.GridLineStyle = ':';  
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
hold on;
 for verline=1:length(timesfunchange)
  plot([timesfunchange(verline) timesfunchange(verline)],get(ax,'YLim'),'--m','LineWidth',linethick);
 end
 legend({'Vdot(t)',legstr},'FontSize',16)
%%  Creating a toroid plot
% u=0:0.07:2*pi;
% v=0:0.1:2*pi;
% a=2;
% b=0.5;
% figure;
%  subplot(1,2,2)
% for i=1:size(u,2)-1
%     for j=1:size(v,2)-1
%         plot3([(a+b*cos(v(j)))*cos(u(i)) (a+b*cos(v(j+1)))*cos(u(i))],[(a+b*cos(v(j)))*sin(u(i)) (a+b*cos(v(j+1)))*sin(u(i))],[b*sin(v(j)) b*sin(v(j+1))]);
%         hold on
%     end
% end
% axis equal
%% Stick simulation and phase plot
% 
%  for i=1:(n-1)
%      subplot(2,2,1)
%      plot3([0 x1(i)],[0 y1(i)],[0 z1(i)],'b',[x1(i) x2(i)],[y1(i) y2(i)],[z1(i) z2(i)],'r','LineWidth',2);
%      axis([-2*l1 2*l1 -2*l1 2*l1 -2*l1 2*l1]);
%      
%      subplot(1,2,2)
%      plot3([(a+b*cos(phi(i)))*cos(theta(i)) (a+b*cos(phi(i)))*cos(theta(i))+0.0001],[(a+b*cos(phi(i)))*sin(theta(i)) (a+b*cos(phi(i)))*sin(theta(i))+0.0001],[b*sin(phi(i)) b*sin(phi(i))+0.0001],'r*');     
% %      
%      subplot(2,2,3)
%      plot([time(i) time(i+1)],[effort(i) effort(i+1)],'b');
%      axis([0 0.05 min([0 min(effort)]) max(effort)]);
%      xlabel('Time');
%      ylabel('Torque applied');
%      title('Control Effort');
%      hold on
%      
%      pause(wait);
%  end
 
% % Stick Simulation only
% figure;
% for i=1:n
%      plot3([0 x1(i)],[0 y1(i)],[0 z1(i)],'b',[x1(i) x2(i)],[y1(i) y2(i)],[z1(i) z2(i)],'r','LineWidth',2);
%      axis([-2*l1 2*l1 -2*l1 2*l1 -2*l1 2*l1]);
%      pause(wait/60);
% end    