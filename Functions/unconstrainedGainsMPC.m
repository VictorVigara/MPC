function [L_x, L_w, L_R, L_u] = unconstrainedGainsMPC(H, phi_x, phi_w, gamma, phi_u, I_0, Q, S)
    % Check that H is positive definite:
    [~, info] = chol(H); % Cholesky factorization
    if info ~= 0
        disp('error, DesignMPCUnconstrainedGains: H is not positive definite')
    end

    % 'CDCECC2011_RobustMPCInnovationForm_finalPublished_2308 Pg7
    Lbar_x = (-H^-1)*gamma'*Q*phi_x; 
    Lbar_w = (-H^-1)*gamma'*Q*phi_w; 
    Lbar_R = (H^-1)*gamma'*Q;        
    Lbar_u = (H^-1)*phi_u'*S*I_0;    

    L_x = I_0'*Lbar_x;
    L_w = I_0'*Lbar_w;
    L_R = I_0'*Lbar_R;
    L_u = I_0'*Lbar_u;
end