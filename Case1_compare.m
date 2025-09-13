%% 
clear
close all
clc
rng(42); 
% load('Data/UAV_NominalParameters.mat');
addpath Functions/
addpath Functions/Visulization/
drone = InitializeDrone_Case1(2);
save_file_path = 'Data/';
figure_disturbances_flag = 1;
figure_states_flag = 1;
figure_inputs_flag = 1;
figure_path_flag = 1;
file_save_flag = 0;
use_directfit = 1;
compute_W_S_flag = 1;
kVal = 1.01e-6;         %Units N * s^2 = kg * m
mVal = 1;               %Units kg
omegaMax2Val = 7.79e6;  %Units (rad / s)^2
gVal = 9.81;            %Units kg * m / s^2
UAVparamValues = [kVal mVal gVal omegaMax2Val];
Case1_compare_result = [];
F_wind = [0.0 0.0 0.0 ]; 
% % Constraints
u_min = [-30/57.3; -30/57.3; 0];
u_max = [30/57.3; 30/57.3; 1];
tstep = 0.001;   % Step size for simulation. Has to be less than minimum allowed Ts
warning( 'on' );
odeOpts = odeset( 'RelTol', 1e-3, 'MaxStep', 0.001 ); % ODE Options
u_eq = drone.u_eq; 
Ts = drone.Ts;  % 50ms
drone.start.x = 0.0;
drone.start.y = 0.0;
drone.start.z = -3.0;  % note that z axis aims to ground
xbound = 10;
ybound = 10;
zbound = 10;
% (x - x_target) constraint
drone.F  = [1/xbound 0 0 0 0 0;...
            -1/xbound   0 0 0 0 0;...
            0     1/ybound 0 0 0 0;...
            0    -1/ybound  0 0 0 0;...
            0     0    1/zbound 0 0 0;...
            0     0   -1/zbound 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0;...
            0 0 0 0 0 0
            ];
drone.F_bar = [drone.F + drone.G*drone.K drone.G*drone.E]; 
% Ax = 4; Ay = 8; Az = 2;    %  Amplitude
% p = 2; q = 1; r = 2;        % Frequency ratio

Ax = 4; Ay = 6; Az =8;     % Amptitude
p = 2; q = 1; r = 2;        % Frequency ratio

phi_x = 0; phi_y = pi/2; phi_z = pi/4;  % Phase
w = 0.4;                       % Fundamental frequency
% Calculation of the common cycle
T_x = 2*pi/(p*w);
T_y = 2*pi/(q*w);
T_z = 2*pi/(r*w);
N = lcm(lcm(p,q),r);
T = 1* 2*pi/w *N;
% T = 5;
tvec = 0:Ts:T;
x_ref = zeros(6,length(tvec));

for k = 1:length(tvec)
    t = tvec(k);
    
    % position
    xk = Ax*sin(p*w*t + phi_x);
    yk = Ay*sin(q*w*t + phi_y);
    zk = -3 + Az*sin(r*w*t + phi_z);
    
    % velocity
    vxk = Ax*p*w*cos(p*w*t + phi_x);
    vyk = Ay*q*w*cos(q*w*t + phi_y);
    vzk = Az*r*w*cos(r*w*t + phi_z);
    
    % [x,y,z,vx,vy,vz]
    x_ref(:,k) = [xk; yk; zk; vxk; vyk; vzk];
end

% plot the reference trajectory
% figure; plot3(x_ref(1,:), x_ref(2,:), x_ref(3,:), 'LineWidth',1.5);
% xlabel('x'); ylabel('y'); zlabel('z');
% grid on; axis equal; title('refs');

%% LQR fly
x_now = x_ref(:,1);
x = x_now';
x_ini = x_now - x_ref(:,1);
time    = 0;
states  = x_now';
inputs  = [];
x_list = x_ini';
allContTime = []; % Variable to hold optimization time
disterbances = []; % Variable to hold the measured disterbances
Fwind = [];
hw = waitbar( 0,'Trying to fly...' ); % create waiting bar

for t = 0:Ts:T-Ts
    waitbar( t/T, hw, 'fly...');
    tic;
    F_wind = [0.05 + 0.02*max(min(randn,1),-1) , ...
            2 + 0.5*max(min(randn,1),-1), ...
            0.04*sin(10*pi*t) + 0.02*max(min(randn,1),-1)];
    k = floor(t/Ts)+1;
    x_target = x_ref(:,k);
    u_k = drone.K*(x_now - x_target);
    x_k = x_now - x_target;
    u = u_k + drone.u_eq;
    u(u < u_min) = u_min(u < u_min);
    u(u > u_max) = u_max(u > u_max);
    contTime = toc;
    % the estimate disturbance
    x_kplus1_estimate = drone.A * x_k + drone.B*(u - drone.u_eq);
    % Setup the simulation model & Simulate
    mod = @(t, state) WindNoiseQuadrotorStateFcnBase( state, u, F_wind, UAVparamValues );
    % Each row in x corresponds to the solution at the value returned in the corresponding row of tt.
    [tt, x] = ode23t( mod, t:tstep:t + Ts, x(end,:), odeOpts );
    % The output is simply all the states (6x1)
    x_now = x(end,:)';
    % Print values in a formatted way
    distur_kplus1 = ((x_now- x_ref(:,k+1)) - x_kplus1_estimate)';
    fprintf(['|time: %.4f  ---u_k = [%.4f %.4f %.4f]         ---distur: [%.4f %.4f %.4f %.4f %.4f %.4f]           \n' ...
            '|state_real: [%.4f %.4f %.4f %.4f %.4f %.4f]   ---state_estimate: [%.4f %.4f %.4f %.4f %.4f %.4f]   \n' ...
            '|uav_kplus1: [%.4f %.4f %.4f %.4f %.4f %.4f]                                                        \n' ...
            '----------------------------------------------------------------------------------------------------\n'],...
        t,[u_k(1,1)*(180/pi) u_k(2,1)*(180/pi) u_k(3,1)],distur_kplus1, x_k ,x_kplus1_estimate,x_now);
    time    = [time;    tt(2:end, end)];
    states  = [states;  x(2:end, :)];
    x_list= [x_list; x_k'];
    inputs  = [inputs;  u'.*ones( size( tt, 1 ) - 1, 1 ) ];
    disterbances = [disterbances; distur_kplus1];
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


figure;
hold on; grid on;
plot3(x_ref(1,:), x_ref(2,:), -x_ref(3,:), 'b-', 'LineWidth', 2);
plot3(LQR.states(:,1), LQR.states(:,2), -LQR.states(:,3), 'r--', 'LineWidth', 1.5);
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m] (down)');
title('3D Trajectory Comparison: Reference vs States');
legend('Reference','States');
axis equal;
view(3);



if compute_W_S_flag
    %% Starting Some precomputation
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('[%s] Fly Finished.\n', datestr(now, 'HH:MM:SS'));
    fprintf('     Starting Some precomputation...\n');
    ddata = LQR.disturbance(10:end,:);  % disturbance data
    predisturb.x = ddata(1:end, 1);  
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
%     [W_A,W_theta,~] = generate_polytope(6,0);
%     [W_A,W_theta,~] = generate_polytope(6,1);
    [W_A,W_theta,~] = generate_polytope(6,2); % this method to get more precise W_fit
    
    if use_directfit
        tic;
        W_true_vertices_T = W_true_vertices';
        for i=1:size(W_A,1)
            [~,idx] = max(W_A(i,:)*W_true_vertices_T);
            W_b_new_i = W_A(i,:)*W_true_vertices_T(:,idx);
            W_b_new(i,1) = W_b_new_i;
        end
        W_true_fit = Polyhedron('A',W_A,'b',W_b_new);
        W_time = toc;
        W_true_fit = W_true_fit.minVRep();
        W_true_fit = W_true_fit.minHRep();
        vwf = W_true_fit.volume;
        W_true = Polyhedron('V',W_true_vertices);
%         W_true = W_true.minVRep();
%         W_true = W_true.minHRep();
        vw = W_true.volume;
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
    tic;
    if 1
    S_samplepoints  = MRPISet(Phi, W_true_fit, 1e-2); 
    convhull_index = convhulln(S_samplepoints);
    S_true_vertices_index = unique(convhull_index(:));
    S_true_vertices             = S_samplepoints(S_true_vertices_index,:);
    
    S_true = Polyhedron(S_true_vertices);
    S_true = minVRep(S_true);
    [S_A,S_theta,~] = generate_polytope(6,1);
    S_true_vertices_T = S_true_vertices';
    for i=1:size(S_A,1)
        [~,idx] = max(S_A(i,:)*S_true_vertices_T);
        S_b_new_i = S_A(i,:)*S_true_vertices_T(:,idx);
        S_b_new(i,1) = S_b_new_i;
    end
    S_true_fit = Polyhedron('A',S_A,'b',S_b_new);
    S_time = toc;
    S_true_fit = S_true_fit.minVRep();
    S_true_fit = S_true_fit.minHRep();
%     PlotW(S_true_fit.V(1,1:3),S_true_fit.V(:,1:3) , S_true_fit.V(:,1:3));
%     PlotW(S_true_fit.V(1,4:6),S_true_fit.V(:,4:6) , S_true_fit.V(:,4:6));
    fprintf('direct S \n');
    tic;
    end

    Case1_compare_result.W_true_fit = W_true_fit;
    Case1_compare_result.S_true_fit = S_true_fit;
    Case1_compare_result.S_con = S_con;
end
