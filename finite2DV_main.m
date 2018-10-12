clc
clear all

%% Accel trajectory
Parms.case = 'accel';
t_total = 60*60*24;
tspan = [0 t_total];
omega = 7.291e-5;
Parms.omega = omega;

options = [];
[t,sol] = ode45(@HCWode,tspan,[0 0 0 0 0 0],options,Parms);

accel = [
          3.*omega.^2.*sol(:,1)+2.*omega.*sol(:,5),... % xdot
         -2.*omega.*sol(:,4),... % ydot
         -omega.^2.*sol(:,3),... % zdot
         ];


state = [sol,accel,t];

FiniteTraj = interp1(t,state,1:1:t(end))';
Parms.FiniteTraj = FiniteTraj;



%% optimization attempt

Knots = 10;
Parms.Knots = Knots;
Parms.TimeTotal = t_total;
Parms.SampleTimes = 1:3600:t_total;

[dvar,fval,exitflag,output] = finite2DV(Parms);


%% plot optimization ouput

DV = dvar;
dt = Parms.TimeTotal/Knots;

rInit = [0 0 0]';
vInit = [0 0 0]';
DVTraj = [];

for i = 1:Knots
    [r,v] = CWHPropagator(rInit,vInit+DV(1:3,i),omega,0:dt-1);
    
    rInit = r(1:3,end);
    vInit = v(1:3,end);
    
    DVTraj = horzcat(DVTraj,[r;v]);
    
end

%% Plot
timeVec = 0:1:length(DVTraj);

ticks = timeVec(1:5*3600:t_total);
labels = ticks/3600;

figure(1)
subplot(2,2,1)
plot(t,sol(:,1)/1000,'-',t,sol(:,2)/1000,'-',t,sol(:,3)/1000,'-')
xlabel('time (secs)')
ylabel('Position (km)')
xticklabels(labels)
xticks(ticks)
legend({'x','y','z'})
title('Reference Position')
axis tight

subplot(2,2,2)
plot(t,sol(:,4),'-',t,sol(:,5),'-',t,sol(:,6),'-')
xlabel('time (secs)')
ylabel('Velocity (m/s)')
xticklabels(labels)
xticks(ticks)
legend({'$\dot{x}$','$\dot{y}$','$\dot{z}$'})
title('Reference Velocity')
axis tight

subplot(2,2,3)
plot(1:length(DVTraj),DVTraj(1,:)/1000,'-',1:length(DVTraj),DVTraj(2,:)/1000,'-',1:length(DVTraj),DVTraj(3,:)/1000,'-')
xlabel('time (secs)')
ylabel('Position (km)')
xticklabels(labels)
xticks(ticks)
legend({'x','y','z'})
title('$\Delta V$ Fitted Position')
axis tight

subplot(2,2,4)
plot(1:length(DVTraj),DVTraj(4,:),'-',1:length(DVTraj),DVTraj(5,:),'-',1:length(DVTraj),DVTraj(6,:),'-')
xlabel('time (secs)')
ylabel('Velocity (m/s)')
xticklabels(labels)
xticks(ticks)
legend({'$\dot{x}$','$\dot{y}$','$\dot{z}$'})
title('$\Delta V$ Fitted Velocity')
axis tight

figure(2)
timeVec = 0:t_total/Knots:t_total;
h = plot(timeVec,DV(1:3,:),'o-');

n = 3;
colors = hsv(n);
set(h, {'color'}, num2cell(colors, 2));

ticks = timeVec;

labels = timeVec/3600;

xticklabels(labels)
xticks(ticks)
xlabel('Time (hours)')
ylabel('km')
title('$\Delta V$ Magnitude')
legend({'$\Delta V_x$','$\Delta V_y$','$\Delta V_z$'})
axis tight


%% Differences

rx = abs(DVTraj(1,end)-FiniteTraj(1,end));
ry = abs(DVTraj(2,end)-FiniteTraj(2,end));
rz = abs(DVTraj(3,end)-FiniteTraj(3,end));
rVx = abs(DVTraj(4,end)-FiniteTraj(4,end));
rVy = abs(DVTraj(5,end)-FiniteTraj(5,end));
rVz = abs(DVTraj(6,end)-FiniteTraj(6,end));
