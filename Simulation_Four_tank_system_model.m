clc; clear all;
%% Simulation scenario
t0 = 0.0;           % [s] Initial time
tf = 20*60;         % [s] Final time
deltaT = 10;        % [s] Sample time 
t = [t0:deltaT:tf]; % [s] sample instants
N = length(t);
m10 = 0.0;          % [g] Liquid mass in tank 1 at time t0
m20 = 0.0;          % [g] Liquid mass in tank 2 at time t0
m30 = 0.0;          % [g] Liquid mass in tank 3 at time t0
m40 = 0.0;          % [g] Liquid mass in tank 4 at time t0
F1 = 300;           % [cm3/s] Flow rate from pump 1
F2 = 300;           % [cm3/s] Flow rate from pump 2
F3 = 0;           % [cm3/s] Flow rate from pump 3
F4 = 0;           % [cm3/s] Flow rate from pump 4
x0 = [m10; m20; m30; m40];
u = [F1; F2];
d = [F3; F4];
step_bin = [1; 1; 1; 1];    % Binary selection for steps for every F
steps = [1.1; 1.25; 1.5];   % Step values for each F
want_step = 0;              % Select if step us wanted
Ts = 4;                     % Sample time for discretize the system 

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
F1_ss = 300;
F2_ss = 300;
F3_ss = 250;
F4_ss = 250;
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
% M = [A B; C D]
% p = eig(A); 
% z = eig(M,N);

% ------------------------------------------------------------------------
% Stochastic simulation 
% ------------------------------------------------------------------------

% Process Noise
Q = [20^2 0;0 40^2];
Lq = chol(Q,'lower');
w = Lq*randn(2,N);

% Measurement Noise
R = eye(4);
Lr = chol(R,'lower');
v = Lr*randn(4,N);

%% Simulate the system
if want_step == 0
    steps = 0;
end

% Solve the system of differential equations
iter_sim = round(tf/deltaT)
x = zeros(1,4);

%Discrete simulation without steps
for kk=1:iter_sim
    y(:,kk) = FourTankSystemSensor(x(kk,:),p); % Get height from mass
    z(:,kk) = FourTankSystemSensor(x(kk,:),p); % Get height from mass
    
    %Define simulation time
    t_ini = (kk-1)*deltaT;
    t_fin = kk*deltaT;

    %Store time-stamps
    time(kk) = t_ini;
    [T,X] = ode15s(@FourTankSystem,[t_ini t_fin],x(kk,:)',[],u + w(:,kk),p,d);
    x(kk+1,:) = X(end,:);
end

time(kk+1) = t_fin;
X = []; T = [];
X(:,:,1) = x;
T(:,:,1) = time';

if want_step == 1
    for i=1:length(steps)
        X(:,:,i) = X(:,:,1);
        T(:,:,i) = T(:,:,1);
    end

    for i=1:length(steps)
        % If steps are applied multiply u and d
        u = [steps(i)*step_bin(1)*F1 steps(i)*step_bin(2)*F2];
        d = [steps(i)*step_bin(3)*F3; steps(i)*step_bin(4)*F4];

        for kk=iter_sim:2*iter_sim
            %Define simulation time
            t_ini = (kk-1)*deltaT;
            t_fin = kk*deltaT;
            %Store time-stamps
            T(kk,:,i) = t_ini;
            [T2,X2] = ode15s(@FourTankSystem,[t_ini t_fin],X(kk,:,i)',[],u,p,d);
            X(kk+1,:,i) = X2(end,:);
        end

        T(kk+1,:,i) = t_fin;
    end
end

% help variables
[nT,nX] = size(X);
a = p(1:4,1)';
A = p(5:8,1)';

% Compute the measured variables
H = zeros(nT,4,length(steps));
for j=1:length(steps)
    for i=1:nT
        H(i,:,j) = FourTankSystemSensor(X(i,:,j),p); % Get height from mass
    end
end

% Compute the flows out of each tank
Qout = zeros(nT,4,length(steps));
for j=1:length(steps)
    for i=1:nT
        Qout(i,:,j) = a.*sqrt(2*g*H(i,:,j));
    end
end

%% PLOTS ------------------------------------------------------------------

% Mass plots
max_m = max(max(max(X)))/1000;

for i=1:length(steps)
    figure(1)

    subplot(2,2,1)
    plot(T(:,:,i)/60,X(:,3,i)/1000)
    xlabel("Time (min)")
    ylabel("Mass (kg)")
    title("Tank 3")
    ylim([0 max_m])
    grid on
    hold on

    subplot(2,2,2)
    plot(T(:,:,i)/60,X(:,4,i)/1000)
    xlabel("Time (min)")
    ylabel("Mass (kg)")
    title("Tank 4")
    ylim([0 max_m])
    grid on
    hold on
    
    subplot(2,2,3)
    plot(T(:,:,i)/60,X(:,1,i)/1000)
    xlabel("Time (min)")
    ylabel("Mass (kg)")
    title("Tank 1")
    ylim([0 max_m])
    grid on
    hold on
    
    subplot(2,2,4)
    plot(T(:,:,i)/60,X(:,2,i)/1000)
    xlabel("Time (min)")
    ylabel("Mass (kg)")
    title("Tank 2")
    ylim([0 max_m])
    grid on

    hold on

    % F_all = [F1; F2; F3; F4];
    % S = zeros(length(H),4,length(F_all));
    % 
    % for i=1:4
    %     for k=1:length(F_all)
    %         for j=1:length(H)
    %             S(j,i,k) = (H(j,i) - ys(i))/abs(F_all(k) - step(k)*F_all(k));
    %         end
    %     end
    % end
    % 
    % 
    
    % 
    % figure(2)
    % subplot(2,2,1)
    % plot(T/60,X(:,3)/1000)
    % xlabel("Time (min)")
    % ylabel("Height (m)")
    % title("Tank 3")
    % ylim([0 max_m])
    % grid on
    % 
    % subplot(2,2,2)
    % plot(T/60,X(:,4)/1000)
    % xlabel("Time (min)")
    % ylabel("Height (m)")
    % title("Tank 4")
    % ylim([0 max_m])
    % grid on
    % 
    % subplot(2,2,3)
    % plot(T/60,X(:,1)/1000)
    % xlabel("Time (min)")
    % ylabel("Height (m)")
    % title("Tank 1")
    % ylim([0 max_m])
    % grid on
    % 
    % subplot(2,2,4)
    % plot(T/60,X(:,2)/1000)
    % xlabel("Time (min)")
    % ylabel("Height (m)")
    % title("Tank 2")
    % ylim([0 max_m])
    % grid on
end

