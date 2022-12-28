function [Q,S,phi_x, phi_w, gamma, phi_u, I_0, H_z, H_S, H] = matricesMPC(A, B, G, C, N, Q_z, S)
% 'CDCECC2011_RobustMPCInnovationForm_finalPublished_2308' Pg7
% Reference weight matrix
Q = kron(eye(N), Q_z); 
% Control input weight matrix
S = kron(eye(N), S); 


nu = size(B,2);
% Phix, phiw and gamma matrices
% 'CDCECC2011_RobustMPCInnovationForm_finalPublished_2308' Pg6

c1 = 1; c2 = 2;
for i=1:N
    % Markov parameters
    Hi = C*A^(i-1)*B;
    k1 = 1; k2 = 2;
    c1_d = c1; c2_d = c2;
    for j=1:N      
        if c1_d<(N*2)
            gamma(c1_d:c2_d,k1:k2) = Hi;
            k1 = k1+2;      k2 = k2+2;
            c1_d = c1_d+2;  c2_d = c2_d+2;  
        end
    end
    phi_x(c1:c2,:) = C*A^(i-1);
    phi_w(c1:c2,:) = C*A^(i-1)*G;
    c1 = c1+2; 
    c2 = c2+2; 
end

phi_u = diag(-1*ones(N-1, 1), -1) + diag(ones(N, 1), 0);
phi_u = kron(phi_u, eye(nu));

I_0 = [eye(nu); zeros((nu*N)-nu, nu)]; 

% setup H:
% 'CDCECC2011_RobustMPCInnovationForm_finalPublished_2308' Pg7
H_z = gamma'*Q*gamma; 
H_S = phi_u'*S*phi_u;
H = H_z + H_S; 


end