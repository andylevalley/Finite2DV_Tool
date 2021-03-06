function [dvar,fval,exitflag,output] = finite2DV(Parameters)
% Parameters Required
%%%%%%%%%%%%%%%%%%%%%
% Parameters.Knots      - 1x1, total number of evenaly spaced knots
% Parameters.TimeTotal  - 1x1, total time of trajectory
% Parameters.omega      - 1x1, mean motion (rad/sec)
% Parameters.FiniteTraj - 6xN, finite trajectory that we are matching

%%


num_dvars = Parameters.Knots + 1;
dvar0 = zeros(3,num_dvars);

Parameters.RefTrajectoryNorm

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
[dvar,fval,exitflag,output] = fminunc(@(dvar)f2DV_objfunc(dvar,Parameters),dvar0,options);


TimeTotal = Parameters.TimeTotal;
Knots = Parameters.Knots;
dt = TimeTotal/Knots;
TimeVec = 0:dt:TimeTotal;

dvar = [dvar; TimeVec];

% options = optimoptions('simulannealbnd','Display','iter');
% 
% [dvar,fval,exitflag,output] = simulannealbnd(@(dvar)f2DV_objfunc(dvar,Parameters),dvar0,[],[],options);

end



