
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
compute_W_S = 1;   % Use predefined W and S or compute from flight data

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
x_now   = x'; 
x_ini = x_now - x_goal';
x_list = [];
time    = 0;
inputs  = [];
states  = x; 
allContTime = []; % Variable to hold MPC computation time
disterbances = []; % Variable to hold the measured disterbances
Fwind = []; % Variable to hold the wind force
Ts = drone.Ts;  % 50ms
hw = waitbar( 0,'Trying to fly...' ); % create waiting bar
%% LQR fly
for t = 0:Ts:T-Ts
    waitbar( t/T, hw, 'fly...');
    tic;
    position_now.y = x_now(2,1);
    F_wind = [1.2 + 0.2*max(min(randn,1),-1), ...
    -0.0 + 0.2*sin(3.8*pi*t) + 0.15*max(min(randn,1),-1), ...
    0.04*sin(10*pi*t) + 0.02*max(min(randn,1),-1)];
    x_target = x_goal';
    x_k = x_now - x_target; 
    u_k = drone.K*(x_k);
    u = u_k + drone.u_eq;
    u(u < u_min) = u_min(u < u_min);
    u(u > u_max) = u_max(u > u_max);
    contTime = toc;
    % the estimate disturbance
    x_mpc = drone.A * x_now + drone.B*(u - drone.u_eq);
    % Setup the simulation model & Simulate
    mod = @(t, state) WindNoiseQuadrotorStateFcnBase( state, u, F_wind, UAVparamValues );
    % Each row in x corresponds to the solution at the value returned in the corresponding row of tt.
    [tt, x] = ode23t( mod, t:tstep:t + Ts, x(end,:), odeOpts );
    % The output is simply all the states (6x1)
    x_now = x(end,:)';
    % Print values in a formatted way
    distur = (x_now - x_mpc)';
    fprintf(['time: %.4f,u = [%.4f %.4f %.4f]---distur: [%.4f %.4f %.4f %.4f %.4f %.4f]\n' ...
    ' state_real: [%.4f %.4f %.4f %.4f %.4f %.4f]---state_estimate: [%.4f %.4f %.4f %.4f %.4f %.4f]\n\n '], t,u',distur, x_now ,x_mpc);
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
LQR.time = time;
LQR.states = states;
LQR.inputs = inputs;
LQR.allContTime = allContTime;
LQR.disturbance = disterbances;
LQR.Fwind = Fwind;
LQR.x_list = x_list;

load Case1_W_S_conservative.mat Conservative  % a predefined conservative W and S, you can also compute it from your own data
W_conservative = Conservative.W;
S_conservative = Conservative.S;


if compute_W_S  % use the LQR flight data to compute W and S
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('[%s] Fly Finished.\n', datestr(now, 'HH:MM:SS'));
    fprintf('     Starting Some precomputation...\n');
    ddata = LQR.disturbance((end-50):end,:);  % disturbance data
    predisturb.x = ddata(1:end, 1);  % 200 * Ts = 10s
    predisturb.y = ddata(1:end, 2);  
    predisturb.z = ddata(1:end, 3);  
    predisturb.u = ddata(1:end, 4);
    predisturb.v = ddata(1:end, 5);  
    predisturb.w = ddata(1:end, 6); 
    % set the vertices
    k = convhulln(ddata);
    hull_vertices = unique(k(:));
    W_true_vertices  = ddata(hull_vertices, :);
    fprintf('============================================================\n\n');
    fprintf(' Fitting \alpha W set - W_true .\n');
    [W_A,W_theta,~] = generate_polytope(6,0); % set the norm vectors of W
    % more precise W can set to: generate_polytope(6,1)
    W_samplepoints_T = ddata';
    for i=1:size(W_A,1)
        [~,idx] = max(W_A(i,:)*W_samplepoints_T);
        W_b_new_i = W_A(i,:)*W_samplepoints_T(:,idx);
        W_b_new(i,1) = W_b_new_i;
    end
    W_true_fit = Polyhedron('A',W_A,'b',W_b_new);
    PlotW(ddata(1:end,1:3), W_conservative.V(:,1:3), W_true_fit.V(:,1:3));
    PlotW(ddata(1:end,4:6), W_conservative.V(:,4:6), W_true_fit.V(:,4:6));
    % compute the safe set S_safe 
    % fprintf('============================================================\n\n');
    % fprintf(' Compute set - S_safe .\n');
    % % the safe constraint of error e (in S_hat): (F+GK)e<=1
    % S_safe = Polyhedron('A',(drone.F + drone.G*drone.K),'b',ones([drone.nc,1]));
    % S_safe = computeVRep(S_safe);
    % Compute Tube S
    fprintf('============================================================\n\n');
    fprintf(' Compute Tube set - S .\n');
    Phi = drone.Phi;
    S_samplepoints  = MRPISet(Phi, W_true_fit, 1e-2); % sampled-based invariant set estimation
    [S_A,S_theta,~] = generate_polytope(6,1);
    S_samplepoints_T = S_samplepoints';
    for i=1:size(S_A,1)
        [~,idx] = max(S_A(i,:)*S_samplepoints_T);
        S_b_new_i = S_A(i,:)*S_samplepoints_T(:,idx);
        S_b_new(i,1) = S_b_new_i;
    end
    S_true_fit = Polyhedron('A',S_A,'b',S_b_new);
    S_true_fit = S_true_fit.minVRep();
    S_true_fit = S_true_fit.minHRep();
    DAMPC.S = S_true_fit;
    DAMPC.W = W_true_fit;
    save('Case1_W_S_DAMPC.mat','DAMPC','-mat');  % save the computed W and S for this flight
end
load Case1_W_S_DAMPC.mat Real  % if compute_W_S = 0, load the last time W and S
drone.S = DAMPC.S;
drone.W = DAMPC.W;
S_list = drone.S;
% Compute h_s
Controller = UQRobustMPC_Case1(drone);
[~, ~, h_s, ~] = Controller.solve([ 0 0 0 0 0 0]');

% compute the feasible region of s_0|k with S
% F_N_real = ComFeasibleRegion_RMPC(drone,h_s,Nu);
% PlotW([0 0 0], F_N_real.V(:,1:3), F_N_real.V(:,1:3));
% save('Case1_F_N_real.mat','F_N_real','-mat');

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
hw = waitbar( 0,'DA-TMPC Trying to fly...' ); % create waiting bar
TrueAndNom = [];
ck =[];
x_list = [];
W_update = [drone.W];
tubeupdatetimelist = [];
for t = 0:Ts:T-Ts
    waitbar( t/T, hw, 'DA-TMPC Trying to fly...');
    % Time step you want to update W and S
    if   t == 1 || t == 1.5 || t == 2.5  || t == 3.5 || t == 4.5 
        fprintf('At time t = %.1f seconds, perform W_update.\n', t);
        tic;
        w_stack = disterbances(end-10:end,:);
        % use points in W_conservative to filter w_stack
        valid_idx = all(W_conservative.A * w_stack' <= W_conservative.b + 1e-9, 1);  
        w_valid   = w_stack(valid_idx, :);
        [W_k ,alpha_k, v_k] = recedingW(W_conservative,w_valid);
        fprintf('at time: %.2f: alpha = %f\n',t, alpha_k);
        W_now = W_k;
        W_update = [W_update; W_now];
        S_now = alpha_k*S_conservative + (1 - alpha_k)*(1 - drone.Phi)^(-1)*v_k;
        tubeupdatetime = toc;
        tubeupdatetimelist = [tubeupdatetimelist,tubeupdatetime];
        drone.W = W_now;
        drone.S = S_now;
        % PlotW(w_stack,W_now.V(:,1:3) , W_conservative.V(:,1:3));
        % PlotW(w_stack,W_now.V(:,4:6) , W_conservative.V(:,4:6));
        % update the controller
        Controller = UQRobustMPC_Case1(drone);
    end 

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
    % use next line to set the same wind for LQR and DATMPC fly
    % F_wind = LQR.Fwind(floor(t/Ts+1),:);
    x_k = (x_now - x_target);
    x_target = x_goal';
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
fprintf('DA-TMPC Fly Finished.\n');
close(hw);
DATMPC.time = time;
DATMPC.states = states;
DATMPC.inputs = inputs;
DATMPC.allContTime = allContTime;
DATMPC.disturbance = disterbances;
DATMPC.Fwind = Fwind;
DATMPC.TrueandNom = TrueAndNom;
DATMPC.ck = ck;
DATMPC.x_list = x_list;
DATMPC.S_list = S_list;


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
    F_wind = DATMPC.Fwind(floor(t/Ts+1),:); % use the same wind setting
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



% Figure DATMPC outcomes
if figure_disturbances_flag
    width = 7.16;
    height = 3.5;
    set(gcf, 'Units', 'inches', 'Position', [1 1 width height]);
    set(gcf, 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 width height], ...
        'PaperSize', [width height]);
    hold on;
    plot(DATMPC.time(5:50:end)', DATMPC.disturbance',LineWidth=1);
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z');
    title('DA-TMPC Disturbance');
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
    sgtitle('DATMPC states');
    zflip_states = DATMPC.states;
    zflip_states(:,3) = -zflip_states(:,3);
    zflip_states(:,6) = -zflip_states(:,6);
    plot(DATMPC.time(1:end)', zflip_states',LineWidth=1);
    plot(DATMPC.time(1:end)', 3.3*ones(size(DATMPC.time(1:end)')),LineWidth=1,Color='red');
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z','p_x^{ub}');
    title('DA-TMPC quadrotor states');
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
    sgtitle('DA-TMPC x_k');
    states = DATMPC.x_list;
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
    plot(DATMPC.time(2:end)', [DATMPC.inputs(:,1)*(180/pi) DATMPC.inputs(:,2)*(180/pi)]', LineWidth=1);
    legend('roll', 'pitch');
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    title('DATMPC Inputs - Roll and Pitch');
    grid on;
    % Subplot 2: Thrust
    subplot(2,1,2);
    plot(DATMPC.time(2:end)', DATMPC.inputs(:,3)', LineWidth=1);
    legend('Thrust');
    xlabel('Time (s)');
    ylabel('Thrust');
    title('DATMPC Inputs - Thrust');
    grid on;
    % Add overall title
    sgtitle('DA-TMPC Inputs');
end

Case1_DATMPC_result.LQR = LQR;
Case1_DATMPC_result.drone  = drone;
Case1_DATMPC_result.DATMPC  = DATMPC;
Case1_DATMPC_result.baseline = baseline;



