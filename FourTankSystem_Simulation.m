clear all; clc; clear figures; close all;

%% Loading simulation scenario and parameters
FourTankSystem_SetupVariables; 

% Define cells to save data to plot it
mass_data = {};
sensor_data = {};
qout_output_data = {};

%% ------------------------------------------------------------------------
% Steady State calculation
% -------------------------------------------------------------------------
% Get index and time step in which first tank achieves steady-state
[X_ss, y_ss, z_ss] = getSteadyState(F1_ss, F2_ss, F3_ss, F4_ss, p, u, d, x0, t0, tf);

%% ------------------------------------------------------------------------
% Compute the solution / Simulate Ode15s
% -------------------------------------------------------------------------

% Solve the system of differential equations
[T,X] = ode15s(@FourTankSystem,[t0 tf],x0,[],u,p,d);
y = FourTankSystemSensor(X,p);

% Save data to plot it
% mass_data = data2Plot(T,X,"ode15s",mass_data);
% sensor_data = data2Plot(T,y,"ode15s",sensor_data);

%% ------------------------------------------------------------------------
% Discret-Time Simulation
%-------------------------------------------------------------------------

% Simulate the system
[X_discret, T_discret, y_discret, z_discret] = discreteSimulation(tf, deltaT, u, p, d);

% Save data to plot it
% mass_data = data2Plot(T_discret,X_discret,"Discrete",mass_data);
% sensor_data = data2Plot(T_discret,y_discret,"Discrete",sensor_data);

%% ------------------------------------------------------------------------
% Stochastic simulation
% -------------------------------------------------------------------------

% Simulate the system
[X_stoch, T_stoch, y_stoch, z_stoch] = stochasticSimulation(tf, deltaT, u, p, d, v, w,...
                                       want_step, steps, step_bin);

height_stoch = FourTankSystemSensor(X_stoch,p); % Get tank height from mass output
qout_stoch = height2flow(height_stoch,p);       % Get out flow from height output

% Save data to plot it
mass_data = data2Plot(T_stoch,X_stoch,"Stochastic",mass_data);
sensor_data = data2Plot(T_stoch, y_stoch, "Stochastic", sensor_data);
qout_output_data = data2Plot(T_stoch, qout_stoch, "Stochastic", qout_output_data);

%% ------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

plotData(mass_data, "mass", "Mass output");
plotData(sensor_data, "height", "Height sensor");
plotData(qout_output_data,"flow", "Output flow");