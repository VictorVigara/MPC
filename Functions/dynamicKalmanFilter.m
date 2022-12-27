function [x_hat, y_hat, z_hat, P, R_z, x, y, z, time] = dynamicKalmanFilter(A, B, C, D, Cz, E, G, u, d, Q, R, S, w, v, x_0, tf, ts, steps)

    % Preallocate space for matrices for storing data:
    x = []; % define all states as 0
    y = []; % define all measurment outputs as 0
    z = []; % define all controlled outputs as 0
    x_hat = [];
    y_hat = [];
    z_hat = [];
    w_hat = [];
    time = zeros(tf+1,1,length(steps));
    P = [];
    Q_k = [];
    R_z = [];

    % As this kalman filter is dynamic, kalman filter gains are
    % calculated every iteration. 



    for k=1:length(steps)
        % Initialize the Kalman filter:
        x(1, :, k) = x_0';  % Initial condition
        x_hat(1, :, k) = x_0';
        time(1,1,k) = 0;
        P(:,:,1,k) = eye(length(x_0)); 
        Q_k(:,:,1,k) = Q;

        % Stochastic simulation 
        for i=1:tf
            % ---- State-space system ---- 
            x(i+1,:,k) = A*x(i,:,k)' + B*u(i,:,k)' + E*d(i,:,k)' + G*w(:,i); 
            y(i,:,k) = C*x(i,:,k)' + D*u(i,:,k)' + v(:,i); 
            z(i,:,k) = Cz*x(i,:,k)'; 
    
            % ---- Measurement update ----
            % Measurement output
            y_hat(i,:,k) = C*x_hat(i,:,k)'; 
            % Innovation
            e_k = (y(i,:,k) - y_hat(i,:,k))'; 
            % --- 
            R_e = C*P(:,:,i,k)*C' + R; 
            % --- Kalman filter gain for states
            K_fx = P(:,:,i,k)*C'*R_e^-1; 
            % --- Kalman filter gain for process noise
            K_fw = S*R_e^-1; 

            % State
            x_hat_k = x_hat(i,:,k)' + K_fx*e_k; 
            % State covariance
            P(:,:,i,k) = P(:,:,i,k) - K_fx*R_e*K_fx';
            % Process noise
            w_hat(:,i,k) = K_fw*e_k; 
            % Process noise covariance
            Q_k(:,:,i,k) = Q - K_fw*R_e*K_fw'; 
    
            % ---- One step prediction (Time update) ----
            time(i+1,1,k) = time(i,1,k) + ts; 
            % State prediction
            x_hat(i+1,:,k) = A*x_hat_k + B*u(i,:,k)' + E*d(i,:,k)' + G*w_hat(:,i,k); 
            % State covariance
            P(:,:,i+1,k) = A*P(:,:,i,k)*A' + G*Q_k(:,:,i,k)*G' - A*K_fx*S'*G' ...
                           - G*S*K_fx'*A'; 
            % Output prediction
            z_hat(i+1,:,k) = Cz*x_hat(i+1,:,k)'; 
            % Output covariance
            R_z(:,:,i+1,k) = Cz*P(:,:,i+1,k)*Cz'; 
        end
    end
    
end