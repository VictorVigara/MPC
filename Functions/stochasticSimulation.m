function [X, T, y, z, u_stoch, d_stoch] = stochasticSimulation(tf, deltaT, u, p, d, v, w, ...
                                                 want_step, steps, step_bin, tfinsteps)    
    % ---- Solve the system of differential equations ----
    iter_sim = round(tf/deltaT);
    Xk = zeros(1,4);
    Tk = 0;
    yk = zeros(1,4);
    zk = zeros(1,4);
    u_stoch(1,:) = u'; 
    d_stoch(1,:) = d'; 

    % If discrete simulation is wanted, adecuate v and w inputs
    if v == 0
        v = zeros(4,tfinsteps);
    end

    if w == 0
        w = zeros(2,tfinsteps);
    end

    % If not steps wanted, set steps to zero
    if want_step == 0
        steps = 0;
    end

    % Loop discrete simulation
    for i=1:iter_sim
        %Define simulation Tk
        t_ini = (i-1)*deltaT;
        t_fin = i*deltaT;

        % Define input with noise
        u_stoch(i+1,:) = (u + w(:,i))'; 
        d_stoch(i+1,:) = d'; 
    
        %Store Tk-stamps
        [~,xk] = ode15s(@FourTankSystem,[t_ini t_fin],Xk(i,:)',[],u_stoch(i+1,:)',p,d_stoch(i+1,:)');
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
            u_stoch(:,:,i) = u_stoch(:,:,1); 
            d_stoch(:,:,i) = d_stoch(:,:,1); 
        end


        for i=1:length(steps)
            % If steps are applied multiply u and d
            if step_bin(1) == 0
                u_step(1,1) = u(1);
            else 
                u_step(1,1) = steps(i)*step_bin(1)*u(1);
            end

            if step_bin(2) == 0
                u_step(1,2) = u(2);
            else
                u_step(1,2) = steps(i)*step_bin(2)*u(2);
            end

            if step_bin(3) == 0
                d_step(1,1) = d(1);
            else
                d_step(1,1) = steps(i)*step_bin(3)*d(1);
            end

            if step_bin(4) == 0
                d_step(2,1) = d(2);
            else
                d_step(2,1) = steps(i)*step_bin(4)*d(2);
            end

%             u_step = [steps(i)*step_bin(1)*u(1) steps(i)*step_bin(2)*u(2)];
%             d_step = [steps(i)*step_bin(3)*d(1); steps(i)*step_bin(4)*d(2)];
            
            for j=iter_sim:tfinsteps
                %Define simulation time
                t_ini = (j)*deltaT;
                t_fin = (j+1)*deltaT;
                %Store time-stamps

                u_stoch(j+2,:,i) = (u_step + w(:,j)'); 
                d_stoch(j+2,:,i) = d_step';
                try
                    [T2,X2] = ode15s(@FourTankSystem,[t_ini t_fin],X(j,:,i)',[],u_stoch(j+2,:,i)',p,d_stoch(j+2,:,i)');
                catch
                    fprintf("Error")
                end
                X(j+2,:,i) = X2(end,:);
                T(j+2,:,i) = t_ini + deltaT;
                y(j+2,:,i) = FourTankSystemSensor(X(j+1,:,i),p) + v(:,j)';  % Get height from mass
                z(j+2,:,i) = FourTankSystemSensor(X(j+1,:,i),p);            % Get height from mass
            end
        end
    end

end