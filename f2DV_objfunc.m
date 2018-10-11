function fval = f2DV_objfunc(dvar,Parms)

Knots = Parms.Knots;
TimeTotal = Parms.TimeTotal;
omega = Parms.omega;


% [x;y;z;xdot;ydot;zdot] (6xN)
FiniteTraj = Parms.FiniteTraj;

rInit = [0 0 0]';
vInit = [0 0 0]';
DVTraj = [];

DV = dvar;


dt = TimeTotal/Knots;
for i = 1:Knots
    [r,v] = CWHPropagator(rInit,vInit+DV(1:3,i),omega,0:dt-1);
    
    rInit = r(1:3,end);
    vInit = v(1:3,end);
    
    DVTraj = horzcat(DVTraj,[r;v]);
    
end

DVTraj(4:6,end) = vInit + DV(1:3,Knots+1);

fval = sum((DVTraj(1,:)-FiniteTraj(1,:)).^2 + (DVTraj(2,:)-FiniteTraj(2,:)).^2 + (DVTraj(3,:)-FiniteTraj(3,:)).^2 +...
           (DVTraj(4,:)-FiniteTraj(4,:)).^2 + (DVTraj(5,:)-FiniteTraj(5,:)).^2 + (DVTraj(6,:)-FiniteTraj(6,:)).^2);

end