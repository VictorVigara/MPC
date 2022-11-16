clear all; clc; clear figures; close all;

%% Loading simulation scenario and parameters
FourTankSystem_SetupVariables; 

% Define cells to save data to plot it
% -Output will be called to variables from differential equation solutions
% -Data will be called to variables from measurements
% -If steps wanted to be saved, step will be used in the variable name
mass_discrete_steps = {};
mass_discrete = {};

%% ------------------------------------------------------------------------
% Steady State calculation
% -------------------------------------------------------------------------
% Get index and time step in which first tank achieves steady-state
[X_ss, y_ss, z_ss] = getSteadyState(F1_ss, F2_ss, F3_ss, F4_ss, p, u, d, x0, t0, Tf);

%% ------------------------------------------------------------------------
% Linearization 
% -------------------------------------------------------------------------
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

%% ------------------------------------------------------------------------
% Discretization of a linear model
% -------------------------------------------------------------------------
[Abar, Bbar] = c2dzoh(A, B, Ts);

%% ------------------------------------------------------------------------
% Transfer function of the linear model
% -------------------------------------------------------------------------
% Obtain transfer functions from linearized differential equations
% [G11_lin ,G12_lin ,G21_lin ,G22_lin] = LDE2tf(A,B,C); 
% M = [A B; C zeros(2,2)];
% N = [];
% poles_lin = eig(A);
% zeros_lin = eig(M,N);

%% ------------------------------------------------------------------------
% Compute the solution / Simulate Ode15s
% -------------------------------------------------------------------------

% Solve the system of differential equations
[T,X] = ode15s(@FourTankSystem,[t0 Tf],x0,[],u,p,d);
y = FourTankSystemSensor(X,p);

% Save data to plot it


%% ------------------------------------------------------------------------
% Discret-Time Simulation
%-------------------------------------------------------------------------

% Simulate stochastic system without error -> v = 0 and w = 0 and without
% steps
[X_discret, T_discret, y_discret, z_discret] = ...
    stochasticSimulation(Tf, deltaT, u, p, d, 0, 0, 0, steps, step_bin);

% Simulate stochastic system without error with steps -> v = 0 and w = 0
[X_discret_steps, T_discret_steps, y_discret_steps, z_discret_steps] = ...
    stochasticSimulation(Tf, deltaT, u, p, d, 0, 0, 1, steps, step_bin);

height_discret = FourTankSystemSensor(X_discret,p); % Get tank height from mass output
qout_discret = height2flow(height_discret,p);       % Get out flow from height output

% Save data to plot it
mass_discrete_steps = data2Plot(T_discret_steps,X_discret_steps,"Deterministic model - ",mass_discrete_steps);
mass_discrete = data2Plot(T_discret,X_discret, "Deterministic model",mass_discrete);

%% ------------------------------------------------------------------------
% Stochastic simulation
% -------------------------------------------------------------------------

% Simulate the system
[X_stoch, T_stoch, y_stoch, z_stoch] = stochasticSimulation(Tf, deltaT, u, p, d, v, w,...
                                       0, steps, step_bin);

height_stoch = FourTankSystemSensor(X_stoch,p); % Get tank height from mass output
qout_stoch = height2flow(height_stoch,p);       % Get out flow from height output

% Save data to plot it
mass_discrete = data2Plot(T_stoch,X_stoch, "Stochastic model",mass_discrete);


%% ------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

plotSteps(mass_discrete_steps, "mass", steps);
plotData(mass_discrete, "mass","Deterministic simulation")