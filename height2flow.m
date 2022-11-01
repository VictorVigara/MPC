function qout = height2flow(h,p)
    % Unpack parameters
    a = p(1:4,1)';
    g = p(11,1);
    
    % Calculate output flow rate for each tanks
    qout = a.*sqrt(2*g.*h);
end