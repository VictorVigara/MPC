function [X, T, y, z] = stochasticSimulation(tf, deltaT, u, p, d, v, w, ...
                                                 want_step, steps, step_bin)    
    % Solve the system of differential equations
    iter_sim = round(tf/deltaT)
    Xk = zeros(1,4);
    Tk = 0;
    yk = zeros(1,4);
    zk = zeros(1,4);

    % Loop discrete simulation
    for i=1:iter_sim
        %Define simulation Tk
        t_ini = (i-1)*deltaT;
        t_fin = i*deltaT;
    
        %Store Tk-stamps
        [~,xk] = ode15s(@FourTankSystem,[t_ini t_fin],Xk(i,:)',[],u + w(:,i),p,d);
        Xk(i+1,:) = xk(end,:);
        Tk(i+1,:) = t_ini + deltaT;
        yk(i+1,:) = FourTankSystemSensor(Xk(i+1,:),p) + v(:,i)';    % Get height from sensor
        zk(i+1,:) = FourTankSystemSensor(Xk(i+1,:),p);              % Get height from mass
    end

    if want_step == 0
        X = Xk;
        T = Tk;
        y = yk;
        z = zk;
    end

    if want_step == 1
        % Initialize variables to store step values
        X = zeros(size(Xk,1), size(Xk,2), length(steps));
        T = zeros(size(Tk,1), size(Tk,2), length(steps));
        y = zeros(size(yk,1), size(yk,2), length(steps));
        z = zeros(size(zk,1), size(zk,2), length(steps));
        
        % Fill output vectors with data before step
        for i=1:length(steps)
            X(:,:,i) = Xk(:,:,1);
            T(:,:,i) = Tk(:,:,1);
            y(:,:,i) = yk(:,:,1);
            z(:,:,i) = zk(:,:,1);
        end

        for i=1:length(steps)
            % If steps are applied multiply u and d
            u_step = [steps(i)*step_bin(1)*u(1) steps(i)*step_bin(2)*u(2)];
            d_step = [steps(i)*step_bin(3)*d(1); steps(i)*step_bin(4)*d(2)];
            
            for j=iter_sim:2*iter_sim
                %Define simulation time
                t_ini = (j)*deltaT;
                t_fin = (j+1)*deltaT;
                %Store time-stamps
                [T2,X2] = ode15s(@FourTankSystem,[t_ini t_fin],X(j,:,i)',[],u_step + w(:,j),p,d_step);
                X(j+2,:,i) = X2(end,:);
                T(j+2,:,i) = t_ini + deltaT;
                y(j+2,:,i) = FourTankSystemSensor(X(j+1,:,i),p) + v(:,j)';  % Get height from mass
                z(j+2,:,i) = FourTankSystemSensor(X(j+1,:,i),p);            % Get height from mass
            end
        end
    end

end