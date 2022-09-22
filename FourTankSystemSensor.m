function ys = FourTankSystemSensor(xs,p)
    % Unpack states, MVs, and parameters
    m = xs; 
    A = p(5:8,1)'; % Tank cross sectional areas [cm2]
    rho = p(12,1); % Density of water [g/cm3]
    ys = m./(rho*A);
end