clc
clear all

%% Accel trajectory
Parms.case = 'accel';
t_total = 60*60*24*1;
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

Knots = 5;
Parms.Knots = Knots;
Parms.TimeTotal = t_total;


% num_dvars = Knots*3;
% dvar0 = zeros(1,num_dvars);

% [dvar,fval,exitflag,output] = fminunc(@(dvar)f2DV_objfunc(dvar,Parms),dvar0);

[dvar,fval,exitflag,output] = finite2DV(Parms);


%% plot optimization ouput

DV = reshape(dvar,[Knots,3])';
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

t = 1:length(Parms.FiniteTraj);
sol = Parms.FiniteTraj;

figure(1)
subplot(2,2,1)
a = plot(t,sol(1,:)/1000,'-',t,sol(2,:)/1000,'-',t,sol(3,:)/1000,'-');
title('Reference Trajectory Position')
xlabel('Time (secs)')
ylabel('Position (km)')
xticklabels(labels)
xticks(ticks)
legend({'x','y','z'})
axis tight

subplot(2,2,2)
b = plot(t,sol(4,:),'-',t,sol(5,:),'-',t,sol(6,:),'-');
title('Reference Trajectory Velocity')
xlabel('Time (secs)')
ylabel('Velocity (m/s)')
xticklabels(labels)
xticks(ticks)
legend({'x','y','z'})
axis tight

subplot(2,2,3)
c = plot(1:length(DVTraj),DVTraj(1,:)/1000,'-',1:length(DVTraj),DVTraj(2,:)/1000,'-',1:length(DVTraj),DVTraj(3,:)/1000,'-');
title('$\Delta$V Fitted Trajectory Position')
xlabel('Time (secs)')
ylabel('Position (km)')
xticklabels(labels)
xticks(ticks)
legend({'x','y','z'})
axis tight

subplot(2,2,4)
d = plot(1:length(DVTraj),DVTraj(4,:),'-',1:length(DVTraj),DVTraj(5,:),'-',1:length(DVTraj),DVTraj(6,:),'-');
title('$\Delta$V Fitted Trajectory Velocity')
xlabel('Time (secs)')
ylabel('Velocity (m/s)')
xticklabels(labels)
xticks(ticks)
legend({'x','y','z'})
axis tight



figure(2)
timeVec = 0:t_total/Knots:t_total;
h = plot(timeVec(1:Knots),DV,'o-');
legend({'$\Delta V_x$','$\Delta V_y$','$\Delta V_z$'});
n = 3;
colors = hsv(n);
set(h, {'color'}, num2cell(colors, 2));
ticks = timeVec;
labels = timeVec/3600;
title('$\Delta$V Magnitude')
xticklabels(labels)
xticks(ticks)
xlabel('Time (hours)')
ylabel('km')
axis([0 t_total -.2 2])

%% End point value

rx = abs(DVTraj(1,end) - sol(1,end));
ry = abs(DVTraj(2,end) - sol(2,end));
rz = abs(DVTraj(3,end) - sol(3,end));
rvx = abs(DVTraj(4,end) - sol(4,end));
rvy = abs(DVTraj(5,end) - sol(5,end));
rvz = abs(DVTraj(6,end) - sol(6,end));

