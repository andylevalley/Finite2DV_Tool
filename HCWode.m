function dxdt = HCWode(t,x,Parms)
persistent n

if isempty(n) == 1
    n = 1;
end

omega = Parms.omega;

        
if t < 1;
    a = [1 1 0];
else
    a = [0 0 0];
end

dxdt = [
        x(4) % x
        x(5) % y
        x(6) % z
        3*omega*omega*x(1)+2*omega*x(5) + a(1) % xdot
        -2*omega*x(4) + a(2) % ydot
        -omega^2*x(3)  + a(3) % zdot
        ];
    
end