function [X, T, y, z, u_rel, d_rel] = linearizedSimulation(Ad, Bd, Cd, Czd, Ed, u_ss, d_ss, tf, deltaT, u,...
                                             p, d, v, w, want_step, steps, step_bin, tfinsteps)
    % Solve the system of differential equations
    iter_sim = round(tf/deltaT);
    Xk = zeros(1,4);
    Tk = 0;
    yk = zeros(1,4);
    zk = zeros(1,2);

    u_relative = u-u_ss;
    d_relative = d-d_ss;

    u_rel(1,:) = u_relative'; 
    d_rel(1,:) = d_relative';

    % If discrete simulation is wanted, adecuate v and w inputs
    if v == 0
        v = zeros(4,2*iter_sim);
    end

    if w == 0
        w = zeros(2,2*iter_sim);
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

        % Define input 
        u_rel(i+1,:) = u_relative'; 
        d_rel(i+1,:) = d_relative';

        % Linear discretized simulation 
        Xk(i+1,:) = (Ad*Xk(end,:)')' + (Bd*u_relative)' + (Ed*d_relative)';
        yk(i+1,:) = (Cd*Xk(end,:)')';   
        zk(i+1,:) = (Czd*Xk(end,:)')';

        % Save simulations results
        Tk(i+1,:) = t_ini + deltaT;
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
            % Save u and d used for simulation
            u_rel(:,:,i) = u_rel(:,:,1);
            d_rel(:,:,i) = d_rel(:,:,1);
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
            
            u_step_rel = u_step'-u;
            d_step_rel = d_step-d;
            for j=iter_sim:tfinsteps
                %Define simulation time
                t_ini = (j)*deltaT;
                t_fin = (j+1)*deltaT;

                % Linear discretized simulation 
                X(j+2,:,i) = (Ad*X(j+1,:,i)' + Bd*u_step_rel + Ed*d_step_rel)';
                y(j+2,:,i) = (Cd*X(j+2,:,i)')';  
                z(j+2,:,i) = (Czd*X(j+2,:,i)')';   

                T(j+2,:,i) = t_ini + deltaT;
                u_rel(j+2,:,i) = u_step_rel';
                d_rel(j+2,:,i) = d_step_rel';
            end
        end
    end

end