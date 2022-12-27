function [x_hat, y_hat, z_hat, x, y, z, time] = statKalmanFilter(A, B, C, D, Cz, E, G, u, d, Q, R, S, w, v, x_0, tf, ts, steps)

    % Preallocate space for matrices for storing data:
    x = []; % define all states as 0
    y = []; % define all measurment outputs as 0
    z = []; % define all controlled outputs as 0
    x_hat = [];
    y_hat = [];
    z_hat = [];
    w_hat = [];
    time = [];



    % Covariance of state computed as Discrete Algebraic Riccati Equation
    P = dare(A',C',G*Q*G',R,G*S); 

    % As this kalman filter is stationary, kalman filter gains are
    % calculated offline. 

    % --- 
    R_e = C*P*C' + R; 
    % --- Kalman filter gain for states
    K_fx = P*C'*R_e^-1; 
    % --- Kalman filter gain for process noise
    K_fw = S*R_e^-1; 

    for k=1:length(steps)
        % Initialize the Kalman filter:
        x(1, :, k) = x_0';  % Initial condition
        x_hat(1, :, k) = x_0';
        time(1,:,k) = 0;

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
            % State
            x_hat_k = x_hat(i,:,k)' + K_fx*e_k; 
            % Process noise
            w_hat(:,i,k) = K_fw*e_k; 
    
            % ---- One step prediction (Time update) ----
            time(i+1,:,k) = time(i,:,k) + ts; 
            % State prediction
            x_hat(i+1,:,k) = A*x_hat_k + B*u(i,:,k)' + E*d(i,:,k)' + G*w_hat(:,i,k); 
            % Output prediction
            z_hat(i+1,:,k) = Cz*x_hat(i+1,:,k)'; 
        end
    end
    
end