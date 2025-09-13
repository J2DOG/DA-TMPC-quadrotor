classdef UQRobustMPC_Case1 < handle
    properties (SetAccess = public)
        K;
        u_eq;
        Px;
        nx;
        nu;
        nc;
        F;
        G;
        F_bar;
        Psi;
        Pc;
        W;
        W_A;
        W_b;
        S;
        SA;
        Sb;
        num_half_space_S;
        N;
        Nu;
        Com_hs;
        UQ;
        optimizer;
    end
    methods (Access = public)
        function obj  = UQRobustMPC_Case1(drone)
            obj.K = drone.K;
            obj.u_eq = drone.u_eq;
            obj.Px = drone.Px;
            obj.nx = drone.nx;
            obj.nu = drone.nu;
            obj.nc = drone.nc;
            obj.F = drone.F;
            obj.G = drone.G;

            obj.F_bar = drone.F_bar;
            obj.Psi = drone.Psi;
            obj.Pc = drone.Pc;
            obj.W     = drone.W;
            obj.W_A   = drone.W.A;
            obj.W_b   = drone.W.b;
            obj.S     = drone.S;
            obj.num_half_space_S = length(drone.S.b);
            obj.SA   = drone.S.A;
            obj.Sb   = drone.S.b;

            obj.N = drone.N;
            obj.Nu = 20;    %TODO
            obj.Com_hs = obj.LP( );
            %obj.UQ = obj.LPforUQ( );
            obj.optimizer = obj.MPCFormulation( );
            


        end  

        function u = solveLQR(obj, x_k) 
            u = obj.K * x_k ;
        end

        function [s_k_opt, u_k, hk, c_k_opt, solvetime] = solve(obj, x_k)
            % solve the u_k = K*x_k + c_0 at time step k
            F_Com = obj.F + obj.G*obj.K;
            hk    = zeros(obj.nc, 1);
            S_A        = obj.SA;
            S_b        = obj.Sb;
%             size(S_A)
%             size(S_b)
            for i = 1:1:obj.nc
%                 fprintf('compute hs: %.4f\n ', i);
                hk(i) = full(obj.Com_hs(F_Com(i, :), S_A, S_b)); % h_k^*
            end
            % Nu_k = Com_Nuk(obj, hk); % this we have checked is equal to Nu (\nu_s), so not used
            tic;
            [s_k_opt, c_k_opt] = obj.optimizer(x_k, hk, S_A, S_b);
            s_k_opt = full(s_k_opt);
            c_k_opt = full(c_k_opt);
            u_k     = obj.K*x_k + c_k_opt(1:obj.nu);
            solvetime = toc;

        end

        function Com_hs = LP(obj)
            opti  = casadi.Opti( );
            F_com = opti.parameter(1, obj.nx);

            S_A   = opti.parameter(obj.num_half_space_S, obj.nx);
            S_b   = opti.parameter(obj.num_half_space_S, 1);

            e     = opti.variable(obj.nx, 1);
            
            opti.minimize(-F_com*e);
            opti.subject_to(S_A*e <= S_b);
            
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            Com_hs = opti.to_function('f', {F_com, S_A, S_b}, {F_com*e});
        end

        function UQ = LPforUQ(obj) % Linear Programming for Updating the Uncertainty Set Online
            opti = casadi.Opti( );
            beta = opti.variable( );
            y    = opti.variable(2, 1);

            alpha_before = opti.parameter( );
            v_before     = opti.parameter(obj.nx, 1);
            w_new        = opti.parameter(obj.nx, 1);

            opti.minimize(-beta);

            opti.subject_to(0 <= beta <= 1);
            opti.subject_to(obj.W_A*y - beta*obj.W_b <= 0);
            opti.subject_to(alpha_before*obj.W_b + (1-alpha_before)*(obj.W_A*v_before) <= (1 - beta)*obj.W_b + obj.W_A*y);
            opti.subject_to(obj.W_A*w_new <= (1 - beta)*obj.W_b + obj.W_A*y);
            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);
            UQ   = opti.to_function('f', {alpha_before, v_before, w_new}, {beta, y});
        end

        function optimizer = MPCFormulation(obj)
            opti = casadi.Opti( );

            s_0k = opti.variable(obj.nx, 1);
            c_k  = opti.variable(obj.N*obj.nu, 1);

            x_k  = opti.parameter(obj.nx, 1);
            hk   = opti.parameter(obj.nc, 1);
            S_A  = opti.parameter(obj.num_half_space_S, obj.nx);
            S_b  = opti.parameter(obj.num_half_space_S, 1);

            J    = s_0k'*obj.Px*s_0k + c_k'*obj.Pc*c_k;
            
            opti.minimize(J);
            opti.subject_to(S_A*(x_k - s_0k) <= S_b);

            for i = 0:1:obj.Nu   % Note here Nu is sufficient, Nu here is the 'vs' in the paper
                opti.subject_to(obj.F_bar*(obj.Psi^i)*[s_0k; c_k] <= ones(obj.nc, 1) - hk);
            end
            opts    = struct('ipopt',struct('print_level',3),'print_time',false);
            opti.solver('ipopt', opts);

            optimizer = opti.to_function('f', {x_k, hk, S_A, S_b}, {s_0k, c_k});
        end

    end

end


