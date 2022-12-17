function [A, B, C, Cz, D, E] = linearizeFTSM(g, ap, rho, At, X_ss, gamma1, gamma2)

T = sqrt((g*ap.^2*rho)./(2*At.*X_ss')); % Derivative of mdot wrt m

A = [-T(1)    0     T(3)     0  ; ...
       0    -T(2)    0      T(4); ...
       0      0    -T(3)     0  ; ...
       0      0      0     -T(4)];

B = [ rho*gamma1            0          ; ...
      0                 rho*gamma2     ; ...
      0                 rho*(1-gamma2) ; ...
      rho*(1-gamma1)        0          ];

C = diag(1./(At*rho));

Cz = C(1:2,:);

D = zeros(4,2);

E = [0     0;
     0     0; 
     rho   0;
     0      rho];

end