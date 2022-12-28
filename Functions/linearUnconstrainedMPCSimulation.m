function [x, y, z, uhat_k] = linearUnconstrainedMPCSimulation(A, B, C, Cz, E, G, tfinsteps, ...
    x_ss, u_0, Q, R, S, v, w, d ,ref, N, L_x, L_w, L_R, L_u)
   x = []; 
   y = []; 
   z = [];
   xhat = []; 
   yhat = []; 
   zhat = []; 
    
    % Initial state
    x(1,:) = x_ss;
    % Initial control input 
    u_km1 = u_0; 
    % Initial state estimate
    xhat(1,:) = x_ss; 

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

    %METER STEPS!!!!!!!!!!!!
    % AHORA K=1
    k=1;
    % Closed loop simulation
    for i=1:tfinsteps
        % System outputs -- State space model
        y(i,:,k) = C*x(i,:,k)' + v(:,i); 
        z(i,:,k) = Cz*x(i,:,k)'; 
        

        % Anticipated future set-points given the current finite horizon:
        % Reference for measured outputs:
        R_k_temp = ref(:,i:(i+N)-1);
        R_k = R_k_temp(:);

        % 
        [e_k, x_hat_k_k, w_hat_k_k, uhat_k(i,:,k), xhat(i+1,:)] = ...
        unconstrainedStatKalmanFilter(A, B, C, G,xhat(i,:), y(i,:,k), R_k, ...
                                      K_fx, K_fw, L_x, L_w, L_R, L_u, u_km1);
        
        u = uhat_k(i, :,k)';

        % ---- Update system -----
        x(i+1,:,k) = A*x(i,:,k)' + B*u + E*d(i,:,k)'+  G*w(:,i);
        
        
        u_km1 = uhat_k(i, :,k)';

    end




end