function drone = InitializeDrone_Case1(Task)
    addpath('Functions/');
    u_eq = [0; 0; 0.3117];
    drone.Ts = 0.05;
    drone.A = [ 1	0	0	0.05	0	0;
                0	1	0	0	0.05	0;
                0	0	1	0	0	0.05;
                0	0	0	1	0	0;
                0	0	0	0	1	0;
                0	0	0	0	0	1];
    drone.B = [0	-0.0122625000000000	0;
               0.0122625000000000	0	0;
                0	0	-0.0393395000000000;
                0	-0.490500000000000	0;
                0.490500000000000	0	0;
                0	0	-1.57358000000000]; 
    % wind
    if Task == 1
     Q = diag([10 10 10 1 1 1]);
     R = diag([200 200 200]);
    end
    %Lissajous
    if Task == 2
     Q = diag([10 10 10 10 10 10]);
     R = diag([50 50 50]);
    end
     

    [K, drone.Px]  = dlqr(drone.A, drone.B, Q, R);
    K = -K;
    drone.K = K;
    Phi = drone.A + drone.B*drone.K;
    eig(Phi)
    drone.Phi = Phi;
    drone.u_eq = u_eq;
    drone.start.x = 0;
    drone.start.y = 0;
    drone.start.z = 0;
    drone.x_target = [0; 0; 0; 0; 0; 0];

%     drone.W = PreSet.W;
%     drone.W = minHRep(drone.W);
% 
%     drone.S = PreSet.S;
%     drone.S = computeHRep(drone.S);
    phi_bound = (60/180)*pi; % rad
    theta_bound = (60/180)*pi; % rad
    z_bound = 80;
    x_bound = 80;
    y_bound = 80;
    T_eq =  0.3117;


    F = [1/x_bound 0 0 0 0 0;...

         -1/x_bound 0 0 0 0 0;...

         0 1/y_bound 0 0 0 0;...
         0 -1/y_bound 0 0 0 0;...
         0 0 1/z_bound 0 0 0;...
         0 0 -1/z_bound 0 0 0;...
         0 0 0 0 0 0;...
         0 0 0 0 0 0;...
         0 0 0 0 0 0;...
         0 0 0 0 0 0;...
         0 0 0 0 0 0;...
         0 0 0 0 0 0
         ];
    G = [0 0 0;...
         0 0 0;...
         0 0 0;...
         0 0 0;...
         0 0 0;...
         0 0 0;...
         1/phi_bound 0 0;...
         -1/phi_bound 0 0;...
         0 1/theta_bound 0;...
         0 -1/theta_bound 0;...
         0 0 1/(1-T_eq);...
         0 0 -1/(T_eq)];
    nc = size(F,1);

    drone.F = F;
    drone.G = G;
    drone.nc = nc;

    nx = 6;
    nu = 3;
    drone.nx = nx;
    drone.nu = nu;

    N = 20; % the length of steps add c_k
    drone.N = N;

    E     = zeros(nu, nu*N);
    drone.E = E;
    E(1:nu,1:nu)  = eye(nu); 

    Pc = zeros(nu*N, nu*N);
    for i    = 1:1:N
        Pc(((i-1)*nu+1):(i*nu),((i-1)*nu+1):(i*nu)) = drone.B' * drone.Px * drone.B + R;
    end
    drone.Pc  = Pc; 

    M   = diag(ones(N-1, 1), 1); 
    M = kron(M,eye(nu));

    Psi = cell(2, 2);
    Psi{1, 1} = Phi;
    Psi{1, 2} = drone.B*E;
    Psi{2, 1} = zeros(nu*N, nx);
    Psi{2, 2} = M;
    drone.Psi   = cell2mat(Psi);

    drone.F_bar = [F + G*K G*E]; 

    save("Data/drone.mat","drone");
end


