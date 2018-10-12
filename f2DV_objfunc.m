function fval = f2DV_objfunc(dvar,Parms)

Knots = Parms.Knots;
TimeTotal = Parms.TimeTotal;
omega = Parms.omega;


% [x;y;z;xdot;ydot;zdot] (6xN)
FiniteTraj = Parms.FiniteTraj;

rInit = [0 0 0]';
vInit = [0 0 0]';
DVTraj = [];


% DV = reshape(dvar,[Knots+1,3])';
DV = dvar;


dt = TimeTotal/Knots;
for i = 1:Knots
    [r,v] = CWHPropagator(rInit,vInit+DV(1:3,i),omega,0:dt-1);
    
    rInit = r(1:3,end);
    vInit = v(1:3,end);
    
    DVTraj = horzcat(DVTraj,[r;v]);
    
end

DVTraj(4:6,end) = vInit + DV(1:3,Knots+1);

step = 1;

% fval = sum((DVTraj(1,1:step:end)-FiniteTraj(1,1:step:end)).^2 + (DVTraj(2,1:step:end)-FiniteTraj(2,1:step:end)).^2 + (DVTraj(3,1:step:end)-FiniteTraj(3,1:step:end)).^2 +...
%            (DVTraj(4,1:step:end)-FiniteTraj(4,1:step:end)).^2 + (DVTraj(5,1:step:end)-FiniteTraj(5,1:step:end)).^2 + (DVTraj(6,1:step:end)-FiniteTraj(6,1:step:end)).^2);


fval = sum((DVTraj(1,end)/1000-FiniteTraj(1,end)/1000)^2 + (DVTraj(2,end)/1000-FiniteTraj(2,end)/1000)^2 + (DVTraj(3,end)/1000-FiniteTraj(3,end)/1000)^2 +...
           (DVTraj(4,end)-FiniteTraj(4,end))^2 + (DVTraj(5,end)-FiniteTraj(5,end))^2 + (DVTraj(6,end)-FiniteTraj(6,end))^2);
end