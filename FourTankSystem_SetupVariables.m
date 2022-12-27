clc; clear all;
%% ------------------------------------------------------------------------
% Simulation scenario
% -------------------------------------------------------------------------

% Simulation time
t0 = 0.0;           % [s] Initial time
Tf = 20*60;         % [s] Final time
deltaT = 1;        % [s] Sample time 
t = [t0:deltaT:Tf]; % [s] sample instants
N = length(t);
Ts = deltaT;                     % Sample time for discretize the system 
tfinsteps = round(3*N);

% Initial mass
m10 = 0.0;          % [g] Liquid mass in tank 1 at time t0
m20 = 0.0;          % [g] Liquid mass in tank 2 at time t0
m30 = 0.0;          % [g] Liquid mass in tank 3 at time t0
m40 = 0.0;          % [g] Liquid mass in tank 4 at time t0
x0 = [m10; m20; m30; m40];

% Steady-state
ss_threshold = 0.05;

% Flow rates pumps
F1 = 300;           % [cm3/s] Flow rate from pump 1
F2 = 300;           % [cm3/s] Flow rate from pump 2
F3 = 250;             % [cm3/s] Flow rate from pump 3
F4 = 250;             % [cm3/s] Flow rate from pump 4
u = [F1; F2];
d = [F3; F4];

% Pumps flow rate for steady state calculation 
F1_ss = 300;
F2_ss = 300;
F3_ss = 250;
F4_ss = 250;
u_ss = [F1_ss; F2_ss];
d_ss = [F3_ss; F4_ss];

% Set-up steps in each pump flow rate
want_step = 1;              % Select if step us wanted
step_bin = [1; 1; 0; 0];    % Binary selection for steps for every F
steps = [1.1; 1.25; 1.5];   % Step values for each F
if want_step == 0
    steps = 0;
end

%Stochastic simulation
% Process Noise
Q = [20^2 0;0 40^2];
Lq = chol(Q,'lower');
w = Lq*randn(2,tfinsteps);
% Measurement Noise
R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,tfinsteps);
% Correlation between the process and measurement noise:
S_wv = zeros(size(Q, 1), size(R, 1)); 

%% ------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
a1 = 1.2272; %[cm2] Area of outlet pipe 1
a2 = 1.2272; %[cm2] Area of outlet pipe 2
a3 = 1.2272; %[cm2] Area of outlet pipe 3
a4 = 1.2272; %[cm2] Area of outlet pipe 4
ap = [a1; a2; a3; a4]; 
A1 = 380.1327; %[cm2] Cross sectional area of tank 1
A2 = 380.1327; %[cm2] Cross sectional area of tank 2
A3 = 380.1327; %[cm2] Cross sectional area of tank 3
A4 = 380.1327; %[cm2] Cross sectional area of tank 4
At = [A1; A2; A3; A4];
gamma1 = 0.58; % Flow distribution constant. Valve 1
gamma2 = 0.68; % Flow distribution constant. Valve 2
gamm = [gamma1; gamma2];
g = 981; %[cm/s2] The acceleration of gravity
rho = 1.00; %[g/cm3] Density of water
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];
