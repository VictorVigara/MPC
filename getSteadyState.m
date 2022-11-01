function [X_ss, y_ss, z_ss] = getSteadyState(F1_ss, F2_ss, F3_ss, F4_ss, p, u, d, x0, t0, tf)
    xs0 = 5000*ones(4,1);   % Initial guess of steady state
    us = [F1_ss; F2_ss];    % Steady-state inputs
    ds = [F3_ss; F4_ss];
    X_ss = (fsolve(@FourTankSystemWrap,xs0, [], us,p,ds))';   % mass at steady state
    y_ss = FourTankSystemSensor(X_ss,p);                        %height in steady state
    z_ss = FourTankSystemOutput(X_ss,p);                        %height in steady state
    
    fprintf("Steady state\n")
    fprintf("Tank 1: %f cm \nTank2: %f cm \nTank3: %f cm \nTank4: %f cm\n",z_ss(1), z_ss(2), z_ss(3), z_ss(4))
end