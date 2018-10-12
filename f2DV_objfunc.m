function fval = f2DV_objfunc(dvar,Parameters)

Knots = Parameters.Knots;
TimeTotal = Parameters.TimeTotal;
omega = Parameters.omega;


% [x;y;z;xdot;ydot;zdot] (6xN)
FiniteTraj = Parameters.FiniteTraj;

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

SampleTimes = Parameters.SampleTimes;


fval = sum((DVTraj(1,SampleTimes(:))-FiniteTraj(1,SampleTimes(:))).^2 +...
           (DVTraj(2,SampleTimes(:))-FiniteTraj(2,SampleTimes(:))).^2 +...
           (DVTraj(3,SampleTimes(:))-FiniteTraj(3,SampleTimes(:))).^2 +...
           (DVTraj(4,SampleTimes(:))-FiniteTraj(4,SampleTimes(:))).^2 +...
           (DVTraj(5,SampleTimes(:))-FiniteTraj(5,SampleTimes(:))).^2 +...
           (DVTraj(6,SampleTimes(:))-FiniteTraj(6,SampleTimes(:))).^2);
end