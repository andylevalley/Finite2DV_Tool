function fval = f2DV_objfunc(dvar,Parms)

Knots = Parms.Knots;
TimeTotal = Parms.TimeTotal;
omega = Parms.omega;


% [x;y;z;xdot;ydot;zdot] (6xN)
FiniteTraj = Parms.FiniteTraj;

rInit = [0 0 0]';
vInit = [0 0 0]';
DVTraj = [];

DV = reshape(dvar,[Knots,3])';

dt = TimeTotal/(Knots-1);
timeVec = 0:dt:TimeTotal;

for i = 1:Knots-1
    [r,v] = CWHPropagator(rInit,vInit+DV(1:3,i),omega,0:dt);
    
    rInit = r(1:3,end);
    vInit = v(1:3,end);
    
    DVTraj = horzcat(DVTraj,[r;v]);
    
end

[r,v] = CWHPropagator(rInit,vInit+DV(1:3,i+1),omega,0);
rInit = r(1:3,end);
vInit = v(1:3,end);

DVTraj = horzcat(DVTraj(:,1:end-1),[r;v]);




fval = sum((DVTraj(1,:)-FiniteTraj(1,:)).^2 + ...
           (DVTraj(2,:)-FiniteTraj(2,:)).^2 + ...
           (DVTraj(3,:)-FiniteTraj(3,:)).^2 + ...
           (DVTraj(4,:)-FiniteTraj(4,:)).^2 + ...
           (DVTraj(5,:)-FiniteTraj(5,:)).^2 + ...
           (DVTraj(6,:)-FiniteTraj(6,:)).^2);
       


end