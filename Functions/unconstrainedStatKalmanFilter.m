function [e_k, x_hat_k_k, w_hat_k_k, uhat_k, x_hat_kp1_k] = unconstrainedStatKalmanFilter(A, B, C, G,x_hat_k_km1, yk, R_k, K_fx, K_fw, L_x, L_w, L_R, L_u, u_km1)
    % 'CDCECC2011_RobustMPCInnovationForm_finalPublished_2308' Pg7
    % ---- Measurement update ----
    % Measurement output
    y_hat_k = C*x_hat_k_km1'; 
    % Innovation
    e_k = (yk - y_hat_k')'; 
    % State
    x_hat_k_k = x_hat_k_km1' + K_fx*e_k; 
    % Process noise
    w_hat_k_k = K_fw*e_k;

    % ---- Regulator ----
    uhat_k = L_x*x_hat_k_k + L_w*w_hat_k_k + L_R*R_k + L_u*u_km1;

    % ---- Update system ----
    x_hat_kp1_k = A*x_hat_k_k + B*uhat_k + G*w_hat_k_k; 
end