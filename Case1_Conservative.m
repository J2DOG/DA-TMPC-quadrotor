%% 
clear
close all
clc
rng(42); 
% load('Data/UAV_NominalParameters.mat');
addpath Functions/
addpath Functions/Visulization/
drone = InitializeDrone_Case1(1);
save_file_path = 'Data/';
figure_disturbances_flag = 1;
figure_states_flag = 1;
figure_inputs_flag = 1;
figure_path_flag = 0;
file_save_flag = 0;
use_directfit = 1;
compute_W_S = 1;
kVal = 1.01e-6;         %Units N * s^2 = kg * m
mVal = 1;               %Units kg
omegaMax2Val = 7.79e6;  %Units (rad / s)^2
gVal = 9.81;            %Units kg * m / s^2
UAVparamValues = [kVal mVal gVal omegaMax2Val];
Case1_result = [];

F_wind = [0.0 0.0 0.0 ]; 
% % Constraints
u_min = [-60/57.3; -60/57.3; 0];
u_max = [60/57.3; 60/57.3; 1];
T = 5;           % Set the simulation final time
tstep = 0.001;   % Step size for simulation. Has to be less than minimum allowed Ts
warning( 'on' );
odeOpts = odeset( 'RelTol', 1e-3, 'MaxStep', 0.001 ); % ODE Options
% Random start points
drone.start.x = 0.0;
drone.start.y = 0.0;
drone.start.z = -3.0;  % note that z axis aims to ground

x = [ drone.start.x, drone.start.y, drone.start.z, 1, 2, 0.2];
x_goal = [ 5, 5, -3, 0, 0, 0];

% (x - x_goal) constraint
drone.F  = [1/0.5 0 0 0 0 0;...
           -1/6  0 0 0 0 0;...
            0     1/0.5 0 0 0 0;...
            0    -1/6  0 0 0 0;...
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
% refs = [.start.x, drone.start.y, drone.start.z,];
allContTime = []; % Variable to hold optimization time
disterbances = []; % Variable to hold the measured disterbances
Fwind = [];
% Controller = UQRobustMPC_Case1(drone);
Ts = drone.Ts;  % 50ms
hw = waitbar( 0,'Trying to fly...' ); % create waiting bar
%% LQR fly
for t = 0:Ts:T-Ts
    waitbar( t/T, hw, 'fly...');
    tic;
    position_now.y = x_now(2,1);
    F_wind = [3*sin(3*pi*t) + 0.2*max(min(randn,1),-1) , ...
          -3*sin(2*pi*t+pi/3) - 0.2*max(min(randn,1),-1), ...
          0.5*sin(10*pi*t) + 0.1*max(min(randn,1),-1)];
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
%     fprintf('time: %.4f,u = [%.4f %.4f %.4f %.4f], position_truth: [%.4f %.4f %.4f]---speed: [%.4f ]\n ',t ,u', x_now(1:3,1) ,norm(x_now(4:6,1)));
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
%     vwind = [vwind; vel_wind];
    Fwind = [Fwind; F_wind];
end

% fprintf('time: %.4f,u = [%.4f %.4f %.4f %.4f]---distur: [%.4f %.4f %.4f %.4f %.4f %.4f]\n state_real: [%.4f %.4f %.4f %.4f %.4f %.4f]---state_estimate: [%.4f %.4f %.4f %.4f %.4f %.4f]\n--\n ', t,u',distur, x_now ,x_mpc);

close(hw);
LQR.time = time;
LQR.states = states;
LQR.inputs = inputs;
LQR.allContTime = allContTime;
LQR.disturbance = disterbances;
LQR.Fwind = Fwind;
LQR.x_list = x_list;
% PlotDisterbence(LQR.disturbance(:,1:3),[0 0 0]);
% PlotDisterbence(LQR.disturbance(:,4:6),[0 0 0]);
if figure_disturbances_flag
    figure;
    sgtitle('LQR Error');
    plot(LQR.time(5:50:end)', LQR.disturbance',LineWidth=1);
    % legend('roll', 'pitch','Thrust');
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z');
    xlabel('Time (s)');
    ylabel('Error');
    grid on;
end
if figure_states_flag
    figure;
    sgtitle('LQR states');
    zflip_states = LQR.states;
    zflip_states(:,3) = -zflip_states(:,3);
    zflip_states(:,6) = -zflip_states(:,6);
    plot(LQR.time(1:end)', zflip_states',LineWidth=1);
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z');
    xlabel('Time (s)');
    ylabel('states');
    grid on;

    figure;
    sgtitle('LQR x_k');
    states = LQR.x_list;
    plot(0:0.05:0.05*(size(states,1)-1)',states',LineWidth=1);
    legend('p_x', 'p_y','p_z','v_x', 'v_y','v_z');
    xlabel('Time (s)');
    ylabel('states');
    grid on;
%% 
end
if figure_inputs_flag
    figure;
    % Subplot 1: Roll and Pitch angles
    subplot(2,1,1);
    plot(LQR.time(2:end)', [LQR.inputs(:,1)*(180/pi) LQR.inputs(:,2)*(180/pi)]', LineWidth=1);
    legend('roll', 'pitch');
    xlabel('Time (s)');
    ylabel('Angle (degrees)');
    title('LQR Inputs - Roll and Pitch');
    grid on;
    % Subplot 2: Thrust
    subplot(2,1,2);
    plot(LQR.time(2:end)', LQR.inputs(:,3)', LineWidth=1);
    legend('Thrust');
    xlabel('Time (s)');
    ylabel('Thrust');
    title('LQR Inputs - Thrust');
    grid on;
    % Add overall title
    sgtitle('LQR Inputs');
end


%-----------------------------------------------------------------------------------------------
%% Starting Some precomputation
if compute_W_S
fprintf('\n');
fprintf('============================================================\n');
fprintf('[%s] Fly Finished.\n', datestr(now, 'HH:MM:SS'));
fprintf('     Starting Some precomputation...\n');
ddata = LQR.disturbance(10:end,:);  % disturbance data
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
fprintf(' Fitting set - W_true .\n');
[W_A,W_theta,~] = generate_polytope(6,0); % this method to get more precise W_fit 
if use_directfit
    W_true_vertices_T = W_true_vertices';
    for i=1:size(W_A,1)
        [~,idx] = max(W_A(i,:)*W_true_vertices_T);
        W_b_new_i = W_A(i,:)*W_true_vertices_T(:,idx);
        W_b_new(i,1) = W_b_new_i;
    end
    W_true_fit = Polyhedron('A',W_A,'b',W_b_new);
    W_true_fit = W_true_fit.minVRep();
    W_true_fit = W_true_fit.minHRep();
    PlotW(ddata(1:end,1:3), W_true_vertices(:,1:3), W_true_fit.V(:,1:3));
    PlotW(ddata(1:end,4:6), W_true_vertices(:,4:6), W_true_fit.V(:,4:6));
end
fprintf('============================================================\n\n');
fprintf(' Compute set - S_safe .\n');
% the safe constraint of error e (in S_hat) (F+GK)e<=1
S_con = Polyhedron('A',(drone.F + drone.G*drone.K),'b',ones([drone.nc,1]));
S_con = computeVRep(S_con);
% Compute Tube S
fprintf('============================================================\n\n');
fprintf(' Compute Tube set - S .\n');
Phi = drone.Phi;
S_samplepoints  = MRPISet(Phi, W_true_fit, 1e-2); 

convhull_index = convhulln(S_samplepoints);
S_true_vertices_index = unique(convhull_index(:));
S_true_vertices             = S_samplepoints(S_true_vertices_index,:);
S_true = Polyhedron(S_true_vertices);
S_true = minVRep(S_true);
[S_A,S_theta,~] = generate_polytope(6,1);
if use_directfit
    S_true_vertices_T = S_true_vertices';
    for i=1:size(S_A,1)
        [~,idx] = max(S_A(i,:)*S_true_vertices_T);
        S_b_new_i = S_A(i,:)*S_true_vertices_T(:,idx);
        S_b_new(i,1) = S_b_new_i;
    end
    S_true_fit = Polyhedron('A',S_A,'b',S_b_new);
    S_true_fit = S_true_fit.minVRep();
    S_true_fit = S_true_fit.minHRep();
    PlotW(S_samplepoints(:,1:3),S_true_fit.V(:,1:3) , S_true_fit.V(:,1:3));
    PlotW(S_samplepoints(:,4:6),S_true_fit.V(:,4:6) , S_true_fit.V(:,4:6));
end
 Conservative.S = S_true_fit;
 Conservative.W = W_true_fit;
 save('Case1_W_S_conservative.mat','Conservative','-mat');
end
load Case1_W_S_conservative.mat Conservative
drone.S = Conservative.S;
drone.W = Conservative.W;
% Compute h_s
Controller = UQRobustMPC_Case1(drone);
[~, ~, h_s, ~] = Controller.solve([ 0 0 0 0 0 0]');

