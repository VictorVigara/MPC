clear all; clc; clear figures; close all;
set(groot,'defaultLegendInterpreter','latex');
%% Loading simulation scenario and parameters
FourTankSystem_SetupVariables; 
iter_sim = round(Tf/deltaT); % iteration at which simulation without steps ends
% Define cells to save data to plot it
% -Output will be called to variables from differential equation solutions
% -Data will be called to variables from measurements
% -If steps wanted to be saved, step will be used in the variable name

% Discrete plots
mass_discrete_steps = {};
z_discrete_steps_plot = {}; 
y_discrete_steps_normalized_F1_plot = {};
y_discrete_steps_normalized_F2_plot = {};

mass_discrete = {};

% Stochastic plots
mass_stochastic_steps = {};
z_stochastic_steps_plot = {}; 
y_stoch_steps_normalized_F1_plot = {}; 
y_stoch_steps_normalized_F2_plot = {}; 

% Linearized plots
mass_linearized_steps = {};

%% ------------------------------------------------------------------------
% Steady State calculation
% -------------------------------------------------------------------------
% Get index and time step in which first tank achieves steady-state
[X_ss, y_ss, z_ss] = getSteadyState(F1_ss, F2_ss, F3_ss, F4_ss, p, u, d, x0, t0, Tf);

%% ------------------------------------------------------------------------
% Linearization 
% -------------------------------------------------------------------------
[A, B, C, Cz, D, E] = linearizeFTSM(g, ap, rho, At, X_ss, gamma1, gamma2);

%% ------------------------------------------------------------------------
% Discretization of a linear model
% -------------------------------------------------------------------------
[Ad, Bd] = c2dzoh(A, B, Ts);
[Ad, Ed] = c2dzoh(A, E, Ts);
Cd = C;
Czd = Cz;
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

% Simulate stochastic system WITHOUT ERROR -> v = 0 and w = 0 and 
% WITHOUT STEPS
[X_discrete, T_discrete, y_discrete, z_discrete] = ...
    stochasticSimulation(Tf, deltaT, u, p, d, 0, 0, 0, steps, step_bin, tfinsteps);

% Simulate stochastic system WITHOUT ERROR WITH STEPS -> v = 0 and w = 0
[X_discrete_steps, T_discrete_steps, y_discrete_steps, z_discrete_steps] = ...
    stochasticSimulation(Tf, deltaT, u, p, d, 0, 0, 1, steps, step_bin, tfinsteps);

T_discrete_steps_init = T_discrete_steps(iter_sim:end,:,:); % Get time from which steps start
y_discrete_steps_init = y_discrete_steps(iter_sim:end,:,:); % Get y from which steps start

qouT_discrete = height2flow(y_discrete,p);                  % Get out flow from height output

% Get normalized output wrt F1 and F2 and different steps
[y_discrete_steps_normalized_F1, y_discrete_steps_normalized_F2] =...
    normalizedStepResponses(y_discrete_steps_init,y_ss,[F1 F2],[F1_ss F2_ss],steps);

% Save data to plot it
mass_discrete_steps = data2Plot(T_discrete_steps,X_discrete_steps, ...
                                "X(t) - Deterministic model - ",mass_discrete_steps);
mass_discrete = data2Plot(T_discrete,X_discrete,...
                          "X(t) - Deterministic model",mass_discrete);

z_discrete_steps_plot = data2Plot(T_discrete_steps,z_discrete_steps, ...
                                "z(t) - Deterministic model - ",z_discrete_steps_plot);
y_discrete_steps_normalized_F1_plot = data2Plot(T_discrete_steps_init,y_discrete_steps_normalized_F1, ...
                                "Norm Steps wrt F1 - Deterministic model - ",y_discrete_steps_normalized_F1_plot);

y_discrete_steps_normalized_F2_plot = data2Plot(T_discrete_steps_init,y_discrete_steps_normalized_F2, ...
                                "Norm Steps wrt F2 - Deterministic model - ",y_discrete_steps_normalized_F2_plot);

% Get tf from normalised step
[y_ss_F1, time_constant_F1]  = getTF(y_discrete_steps_normalized_F1,steps,T_discrete_steps_init);
evalTF(y_ss_F1, time_constant_F1, steps, y_discrete_steps_normalized_F1, T_discrete_steps_init, "F1");
%% ------------------------------------------------------------------------
% Stochastic simulation
% -------------------------------------------------------------------------

% Simulate the system
[X_stoch, T_stoch, y_stoch, z_stoch] = ...
stochasticSimulation(Tf, deltaT, u, p, d, v, w, 0, steps, step_bin, tfinsteps);

[X_stoch_steps, T_stoch_steps, y_stoch_steps, z_stoch_steps] = ...
stochasticSimulation(Tf, deltaT, u, p, d, v, w,1, steps, step_bin, tfinsteps);

T_stoch_steps_init = T_stoch_steps(iter_sim:end,:,:);
y_stoch_steps_init = y_stoch_steps(iter_sim:end,:,:);

height_stoch = FourTankSystemSensor(X_stoch,p); % Get tank height from mass output
qout_stoch = height2flow(height_stoch,p);       % Get out flow from height output

% Get normalized output wrt F1 and F2 and different steps
[y_stochastic_steps_normalized_F1, y_stochastic_steps_normalized_F2] =...
    normalizedStepResponses(y_stoch_steps_init,y_ss,[F1 F2],[F1_ss F2_ss],steps);

% Save data to plot it
mass_stochastic_steps = data2Plot(T_stoch_steps, X_stoch_steps, "X(t) - Stochastic model",mass_stochastic_steps);

z_stochastic_steps_plot = data2Plot(T_stoch_steps, z_stoch_steps, "z(t) - Stochastic model",z_stochastic_steps_plot);

y_stoch_steps_normalized_F1_plot = data2Plot(T_stoch_steps_init,y_stochastic_steps_normalized_F1, ...
                                "Norm Steps wrt F1 - Stochastic model - ",y_stoch_steps_normalized_F1_plot);

y_stoch_steps_normalized_F2_plot = data2Plot(T_stoch_steps_init,y_stochastic_steps_normalized_F2, ...
                                "Norm Steps wrt F2 - Stochastic model - ",y_stoch_steps_normalized_F2_plot);



%% ------------------------------------------------------------------------
% Linearized simulation
% -------------------------------------------------------------------------
X0 = [0;0;0;0];

[X, T, y, z] = linearizedSimulation(Ad, Bd, Cd, Czd, Ed, u_ss, d_ss, Tf, deltaT, u,...
                                    p, d, v, w, want_step, steps, step_bin, tfinsteps);
X_lin = X+X_ss;
mass_linearized_steps = data2Plot(T, X_lin, 'X(t) - Linearized model', mass_linearized_steps);
%% ------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

plotSteps(mass_discrete_steps, "mass", steps);
plotSteps(mass_stochastic_steps, "mass", steps);
plotSteps(mass_linearized_steps, "mass", steps);
plotSteps(z_stochastic_steps_plot, "height", steps); 
plotSteps(z_discrete_steps_plot, "height", steps);
plotSteps(y_discrete_steps_normalized_F1_plot, "height", steps); 
plotSteps(y_discrete_steps_normalized_F2_plot, "height", steps);
plotSteps(y_stoch_steps_normalized_F1_plot, "height", steps); 
plotSteps(y_stoch_steps_normalized_F2_plot, "height", steps); 

