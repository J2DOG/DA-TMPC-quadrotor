
clear
close all
clc
rng(42); 
addpath Functions/
addpath Functions/Visulization/
drone = InitializeDrone_Case1(1);
save_file_path = 'Data/';

% Flags
figure_disturbances_flag = 1;
figure_states_flag = 1;
figure_inputs_flag = 1;
figure_path_flag = 0;
file_save_flag = 0;
use_directfit = 1; 
compute_W_S = 0;   % Use predefined W and S or compute from flight data

kVal = 1.01e-6;         %Units N * s^2 = kg * m
mVal = 1;               %Units kg
omegaMax2Val = 7.79e6;  %Units (rad / s)^2
gVal = 9.81;            %Units kg * m / s^2
UAVparamValues = [kVal mVal gVal omegaMax2Val];
Case1_result = [];


F_wind = [0.0 0.0 0.0 ]; 
% Input Constraints
u_min = [-60/57.3; -60/57.3; 0];
u_max = [60/57.3; 60/57.3; 1];
T = 5;           % Set the simulation final time
tstep = 0.001;   % Step size for simulation. Has to be less than minimum allowed Ts
warning( 'on' );
odeOpts = odeset( 'RelTol', 1e-3, 'MaxStep', 0.001 ); % ODE Options
% start position of the quadrotor
drone.start.x = 0.0;
drone.start.y = 0.0;
drone.start.z = -3.0;  % note that z axis aims to ground
drone.start.vx = 1.0;
drone.start.vy = 2.0;
drone.start.vz = 0.2;
x = [ drone.start.x, drone.start.y, drone.start.z, drone.start.vx, drone.start.vy, drone.start.vz];
x_goal = [ 5, 5, -3, 0, 0, 0]; % target position of the quadrotor
% (x - x_goal) constraint
drone.F  = [1/0.5 0 0 0 0 0;...
            -1/5  0 0 0 0 0;...
            0     1/0.5 0 0 0 0;...
            0    -1/5  0 0 0 0;...
            0     0    1/5 0 0 0;...
            0     0   -1/5 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0
            ];

drone.F_bar = [drone.F + drone.G*drone.K drone.G*drone.E]; 
drone.nc = size(drone.F,1);
u_eq = drone.u_eq;
Ts= drone.Ts;
load Case1_W_S_conservative.mat Conservative
Conservative.S;
drone.W = Conservative.W;
drone.S = Conservative.S;
%-----------------------------------------------------------------------------------------------
%% TubeMPC Fly
x = [ drone.start.x, drone.start.y, drone.start.z, drone.start.vx, drone.start.vy, drone.start.vz];
x_now   = x';   
time    = 0;
inputs  = [];
states  = x; 
allContTime = []; % Variable to hold optimization time
disterbances = []; % Variable to hold the measured disterbances
vwind =[0 0 0];
Fwind =[0 0 0];
Controller = UQRobustMPC_Case1(drone);
hw = waitbar( 0,'TMPC Trying to fly...' ); % create waiting bar
TrueAndNom = [];
ck =[];
x_list = [];
S_list = [drone.S];
tubeupdatetimelist = [];
for t = 0:Ts:T-Ts
    waitbar( t/T, hw, 'TMPC Trying to fly...');

    % different wind force in different stages   
    if t >= 0 && t < 1
        F_wind = [1.2 + 0.2*max(min(randn,1),-1), ...
                -0.0 + 0.2*sin(3.8*pi*t) + 0.15*max(min(randn,1),-1), ...
                0.04*sin(10*pi*t) + 0.02*max(min(randn,1),-1)];
    elseif t >= 1 && t < 2
        F_wind = [2 + 0.2*max(min(randn,1),-1), ...
                -0.5 + 0.2*sin(3.8*pi*t) + 0.15*max(min(randn,1),-1), ...
                0.04*sin(10*pi*t) + 0.02*max(min(randn,1),-1)];
    elseif t >= 2 && t < 3
        F_wind = [0.9 + 0.2*max(min(randn,1),-1), ...
            -0.2 + 0.2*sin(3.8*pi*t) + 0.15*max(min(randn,1),-1), ...
            0.04*sin(10*pi*t) + 0.02*max(min(randn,1),-1)];
    elseif t >= 3 && t <= (T - Ts)
        F_wind = [1.5 + 0.2*max(min(randn,1),-1), ...
            0.4 + 0.2*sin(3.8*pi*t) + 0.15*max(min(randn,1),-1), ...
            0.04*sin(10*pi*t) + 0.02*max(min(randn,1),-1)];
    else
        F_wind = [0, 0, 0];
    end
    % use next line to set the same wind for LQR and TMPC fly
    % F_wind = LQR.Fwind(floor(t/Ts+1),:);
    x_target = x_goal';
    x_k = (x_now - x_target);
    [s_k_opt, u_k, hk, c_k_opt,contTime] = Controller.solve(x_now - x_target);
    TrueAndNom = [TrueAndNom; x_now',(s_k_opt + x_target)'];
    ck = [ck; c_k_opt(1:3,1)'];
    u = u_k + drone.u_eq;
    % the estimate disturbance
    x_mpc = drone.A * x_k + drone.B*u_k;
    % Setup the simulation model & Simulate
    mod = @(t, state) WindNoiseQuadrotorStateFcnBase( state, u, F_wind, UAVparamValues );
    % Each row in x corresponds to the solution at the value returned in the corresponding row of tt.
    [tt, x] = ode23t( mod, t:tstep:t + Ts, x(end,:), odeOpts );
    % The output is simply all the states (6x1)
    x_now = x(end,:)';
    % Print values in a formatted way
    distur = (x_now-x_target - x_mpc)';
    % Optimized fprintf with better formatting and readability
    fprintf('\n=== Time: %.4f ===\n', t);
    fprintf('Control Input u:     [%.4f, %.4f, %.4f]\n', u');
    fprintf('Optimal c:           [%.4f, %.4f, %.4f]\n', c_k_opt(1:3,1)');
    fprintf('Disturbance:         [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', distur);
    fprintf('Real State:          [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', (x_now - x_target)');
    fprintf('Estimated State:     [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', (x_mpc)');
    fprintf('hk:                  [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', hk');
    fprintf('====================================================================\n\n');
    time    = [time;    tt(2:end, end)];
    states  = [states;  x(2:end, :)];
    inputs  = [inputs;  u'.*ones( size( tt, 1 ) - 1, 1 ) ];
    x_list= [x_list; x_k'];
    disterbances = [disterbances; distur];
    % Save the MPC computation time for this step
    allContTime = [allContTime; contTime];
    Fwind = [Fwind; F_wind];
    S_list =[S_list; drone.S];
end
fprintf('TMPC Fly Finished.\n');
close(hw);
TMPC.time = time;
TMPC.states = states;
TMPC.inputs = inputs;
TMPC.allContTime = allContTime;
TMPC.disturbance = disterbances;
TMPC.Fwind = Fwind;
TMPC.TrueandNom = TrueAndNom;
TMPC.ck = ck;
TMPC.x_list = x_list;
TMPC.S_list = S_list;
% baseline fly (LQR in same wind config)
x = [ drone.start.x, drone.start.y, drone.start.z, 1, 2, 0.2];
x_now   = x';  
time    = 0;
inputs  = [];
states  = x; 
allContTime = []; % Variable to hold optimization time
disterbances = []; % Variable to hold the measured disterbances
Fwind = [];
hw = waitbar( 0,'Trying to fly...' ); % create waiting bar
x_list = [];
for t = 0:Ts:T-Ts
    waitbar( t/T, hw, 'fly...');
    tic;
    position_now.y = x_now(2,1);
    F_wind = TMPC.Fwind(floor(t/Ts+1),:); % use the same wind setting
    x_target = x_goal';
    x_k = x_now - x_target;
    u_k = drone.K*(x_k);
    u = u_k + drone.u_eq;
    u(u < u_min) = u_min(u < u_min);
    u(u > u_max) = u_max(u > u_max);
    contTime = toc;
    % the estimate disturbance
    x_mpc = drone.A * x_k + drone.B*(u - drone.u_eq);
    % Setup the simulation model & Simulate
    mod = @(t, state) WindNoiseQuadrotorStateFcnBase( state, u, F_wind, UAVparamValues );
    % Each row in x corresponds to the solution at the value returned in the corresponding row of tt.
    [tt, x] = ode23t( mod, t:tstep:t + Ts, x(end,:), odeOpts );
    % The output is simply all the states (6x1)
    x_now = x(end,:)';
    distur = (x_now - x_target - x_mpc)';
     time    = [time;    tt(2:end, end)];
    states  = [states;  x(2:end, :)];
    inputs  = [inputs;  u'.*ones( size( tt, 1 ) - 1, 1 ) ];
    x_list= [x_list; x_k'];
    disterbances = [disterbances; distur];
    % Save the computation time for this step
    allContTime = [allContTime; contTime];
    Fwind = [Fwind; F_wind];
end
close(hw);
baseline.time = time;
baseline.states = states;
baseline.inputs = inputs;
baseline.allContTime = allContTime;
baseline.disturbance = disterbances;
baseline.Fwind = Fwind;
baseline.x_list = x_list;



% Figure TMPC outcomes
if figure_disturbances_flag
    width = 7.16;
    height = 3.5;
    set(gcf, 'Units', 'inches', 'Position', [1 1 width height]);
    set(gcf, 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 width height], ...
        'PaperSize', [width height]);
    hold on;
    plot(TMPC.time(5:50:end)', TMPC.disturbance',LineWidth=1);
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z');
    title('TMPC Disturbance');
    xlabel('Time (s)');
    ylabel('Error');
    grid on;
end

if figure_states_flag
    width = 7.16;  
    height = 3.5;
    set(gcf, 'Units', 'inches', 'Position', [1 1 width height]);
    set(gcf, 'PaperUnits', 'inches', ...
            'PaperPosition', [0 0 width height], ...
            'PaperSize', [width height]);
    hold on;
    sgtitle('TMPC states');
    zflip_states = TMPC.states;
    zflip_states(:,3) = -zflip_states(:,3);
    zflip_states(:,6) = -zflip_states(:,6);
    plot(TMPC.time(1:end)', zflip_states',LineWidth=1);
    plot(TMPC.time(1:end)', 3.3*ones(size(TMPC.time(1:end)')),LineWidth=1,Color='red');
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z','p_x^{ub}');
    title('TMPC quadrotor states');
    xlabel('Time (s)');
    ylabel('states');
    hold off;
    grid on;

    width = 7.16;  
    height = 3.5;
    set(gcf, 'Units', 'inches', 'Position', [1 1 width height]);
    set(gcf, 'PaperUnits', 'inches', ...
            'PaperPosition', [0 0 width height], ...
            'PaperSize', [width height]);
    sgtitle('TMPC x_k');
    states = TMPC.x_list;
    plot(0:0.05:0.05*(size(states,1)-1)',states',LineWidth=1);
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z');
    xlabel('Time (s)');
    ylabel('states');
    grid on;
end

if figure_inputs_flag
    figure;
    % Subplot 1: Roll and Pitch angles
    subplot(2,1,1);
    plot(TMPC.time(2:end)', [TMPC.inputs(:,1)*(180/pi) TMPC.inputs(:,2)*(180/pi)]', LineWidth=1);
    legend('roll', 'pitch');
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    title('TMPC Inputs - Roll and Pitch');
    grid on;
    % Subplot 2: Thrust
    subplot(2,1,2);
    plot(TMPC.time(2:end)', TMPC.inputs(:,3)', LineWidth=1);
    legend('Thrust');
    xlabel('Time (s)');
    ylabel('Thrust');
    title('TMPC Inputs - Thrust');
    grid on;
    % Add overall title
    sgtitle('TMPC Inputs');
end
Case1_TMPC_result.drone  = drone;
Case1_TMPC_result.TMPC  = TMPC;
Case1_TMPC_result.baseline = baseline;



