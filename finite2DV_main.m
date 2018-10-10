clc
clear all

%% Accel trajectory
Parms.case = 'accel';
t_total = 60*60*12;
tspan = [0: t_total];
omega = 7.291e-5;
Parms.omega = omega;

options = [];
[t,sol] = ode45(@HCWode,tspan,[0 0 0 0 0 0],options,Parms);

figure(1)
plot(t,sol(:,1)/1000,'-',t,sol(:,2)/1000,'-',t,sol(:,3)/1000,'-')
hold on
plot(t,sol(:,4),'-',t,sol(:,5),'-',t,sol(:,6),'-')

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
Parms.TimeTotal = 60*60*12;


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

figure(2)
plot(1:length(DVTraj),DVTraj(1,:)/1000,'-',1:length(DVTraj),DVTraj(2,:)/1000,'-',1:length(DVTraj),DVTraj(3,:)/1000,'-')
hold on
plot(1:length(DVTraj),DVTraj(4,:),'-',1:length(DVTraj),DVTraj(5,:),'-',1:length(DVTraj),DVTraj(6,:),'-')


