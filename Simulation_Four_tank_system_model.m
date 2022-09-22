%% Simulation scenario
t0 = 0.0; % [s] Initial time
tf = 20*60; % [s] Final time
m10 = 0.0; % [g] Liquid mass in tank 1 at time t0
m20 = 0.0; % [g] Liquid mass in tank 2 at time t0
m30 = 0.0; % [g] Liquid mass in tank 3 at time t0
m40 = 0.0; % [g] Liquid mass in tank 4 at time t0
F1 = 300; % [cm3/s] Flow rate from pump 1
F2 = 300; % [cm3/s] Flow rate from pump 2
x0 = [m10; m20; m30; m40];
u = [F1; F2];
step = 1.25;
want_step = 0;
Ts = 4; % Sample time for discretize the system 

%% Parameters
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
gamma1 = 0.45; % Flow distribution constant. Valve 1
gamma2 = 0.40; % Flow distribution constant. Valve 2
gamm = [gamma1; gamma2];
g = 981; %[cm/s2] The acceleration of gravity
rho = 1.00; %[g/cm3] Density of water
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

%% Get steady-state point
xs0 = 5000*ones(4,1); % Initial guess of steady state
us = [F1; F2]; % Steady-state inputs
xs = fsolve(@FourTankSystemWrap,xs0, [], us,p); % mass at steady state
ys = FourTankSystemSensor(xs',p); %height in steady state
zs = FourTankSystemSensor(xs',p); %height in steady state

% Get time in which steady state is reached
[T,X] = ode15s(@FourTankSystem,[t0 tf],x0,[],u,p);

dif = abs(xs'-X);
mean_err = mean(abs(dif),2);
idx_ss = find(mean_err == (min(mean_err))); % idx for the steady state
t_ss = T(idx_ss);

T =[]; X = [];

%% Linearization
hs = ys;
Ti = (ap./At).*sqrt(g./(2*hs'));
A = [-Ti(1) 0 Ti(3) 0; 0 -Ti(2) 0 Ti(4); 0 0 -Ti(3) 0; 0 0 0 -Ti(4)];
B = [rho*gamm(1) 0; 0 rho*gamm(2); 0 rho*(1-gamm(2)); rho*(1-gamm(1)) 0];
C = diag(1./(rho*At));
Cz = C(1:2,:);

%% Discretization
[Abar, Bbar] = c2dzoh(A, B, Ts);

%% Linear model simulation

%% Pole-Zero representation
M = [A B; C D]
p = eig(A); 
z = eig(M,N);

%% Simulate the system
% Solve the system of differential equations
if want_step == 0
    [T,X] = ode15s(@FourTankSystem,[t0 tf],x0,[],u,p);
else
    [T,X] = ode15s(@FourTankSystem,[t0 t_ss],x0,[],u,p);

    % Step in actuator
    u = [step*F1 F2];
    x0 = X(end,:);
    [T2,X2] = ode15s(@FourTankSystem,[t_ss tf],x0,[],u,p);

    T = [T;T2]; X = [X;X2];
end

% help variables
[nT,nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';
% Compute the measured variables
H = zeros(nT,nX);
for i=1:nT
H(i,:) = FourTankSystemSensor(X(i,:),p); % Get height from mass
end
% Compute the flows out of each tank
Qout = zeros(nT,nX);
for i=1:nT
Qout(i,:) = a.*sqrt(2*g*H(i,:));
end

max_m = max(max(X))/1000;

figure(1)
subplot(2,2,1)
plot(T/60,X(:,3)/1000)
xlabel("Time (min)")
ylabel("Mass (kg)")
title("Tank 3")
ylim([0 max_m])
grid on

subplot(2,2,2)
plot(T/60,X(:,4)/1000)
xlabel("Time (min)")
ylabel("Mass (kg)")
title("Tank 4")
ylim([0 max_m])
grid on

subplot(2,2,3)
plot(T/60,X(:,1)/1000)
xlabel("Time (min)")
ylabel("Mass (kg)")
title("Tank 1")
ylim([0 max_m])
grid on

subplot(2,2,4)
plot(T/60,X(:,2)/1000)
xlabel("Time (min)")
ylabel("Mass (kg)")
title("Tank 2")
ylim([0 max_m])
grid on

