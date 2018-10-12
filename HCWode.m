function dxdt = HCWode(t,x,Parms)

omega = Parms.omega;

        
if t < 20
    a = [.01 .01 .01];
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