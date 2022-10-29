function [Xk, Tk, yk, zk] = stochasticSimulation(tf, deltaT, u, p, d, v, w)    
    % Solve the system of differential equations
    iter_sim = round(tf/deltaT)
    Xk = zeros(1,4);
    Tk = 0;
    yk = zeros(1,4);
    zk = zeros(1,4);
    %Discrete simulation without steps
    for i=1:iter_sim
        
        %Define simulation Tk
        t_ini = (i-1)*deltaT;
        t_fin = i*deltaT;
    
        %Store Tk-stamps
        [~,xk] = ode15s(@FourTankSystem,[t_ini t_fin],Xk(i,:)',[],u + w(:,i),p,d);
        Xk(i+1,:) = xk(end,:);
        Tk(i+1,:) = t_ini + deltaT;
        yk(i+1,:) = FourTankSystemSensor(Xk(i+1,:),p) + v(:,i)'; % Get height from mass
        zk(i+1,:) = FourTankSystemSensor(Xk(i+1,:),p); % Get height from mass
    end