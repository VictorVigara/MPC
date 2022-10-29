clear all; clc; clear figures; close all;
%% Loading simulation scenario and parameters
FourTankSystem_SetupVariables; 

% Define cells to save data to plot it
mass_data = {};
sensor_data = {};

%% ------------------------------------------------------------------------
% Steady State calculation
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Get steady-state point
% -------------------------------------------------------------------------

xs0 = 5000*ones(4,1); % Initial guess of steady state
us = [F1_ss; F2_ss]; % Steady-state inputs
ds = [F3_ss; F4_ss];
xs = fsolve(@FourTankSystemWrap,xs0, [], us,p,ds); % mass at steady state
ys = FourTankSystemSensor(xs',p); %height in steady state
zs = FourTankSystemSensor(xs',p); %height in steady state


% Get time in which steady state is reached
[T,X] = ode15s(@FourTankSystem,[t0 tf],x0,[],u,p,d);

dif = abs(xs'-X);
mean_err = mean(abs(dif),2);
idx_ss = find(mean_err == (min(mean_err))); % idx for the steady state
t_ss = T(idx_ss);

%% ------------------------------------------------------------------------
% Compute the solution / Simulate Ode15s
% -------------------------------------------------------------------------

% Solve the system of differential equations
[T,X] = ode15s(@FourTankSystem,[t0 tf],x0,[],u,p,d);
y = FourTankSystemSensor(X,p);

% Save data to plot it
mass_data = data2Plot(T,X,"ode15s",mass_data);
sensor_data = data2Plot(T,y,"ode15s",sensor_data);

%% ------------------------------------------------------------------------
% Discret-Time Simulation
%-------------------------------------------------------------------------

% Simulate the system
[X_discret, T_discret, y_discret, z_discret] = discreteSimulation(tf, deltaT, u, p, d);

% Save mass data to plot it
mass_data = data2Plot(T_discret,X_discret,"Discrete",mass_data);
% Save sensor data
sensor_data = data2Plot(T_discret,y_discret,"Discrete",sensor_data);

%% ------------------------------------------------------------------------
% Stochastic simulation
% -------------------------------------------------------------------------

% Simulate the system
[X_stoch, T_stoch, y_stoch, z_stoch] = stochasticSimulation(tf, deltaT, u, p, d, v, w)

% Save mass data to plot it
mass_data = data2Plot(T_stoch,X_stoch,"Stochastic",mass_data)
% Save sensor data
sensor_data = data2Plot(T_stoch, y_stoch, "Stochastic", sensor_data)

%% ------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------

plotData(mass_data, "mass", "Mass");
plotData(sensor_data, "height", "Height sensor")
