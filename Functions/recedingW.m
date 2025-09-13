function [W_k, alpha_k, v_k] = recedingW(W,w_stack)
%RECEDINGW Summary of this function goes here
%   W: Polyhedron w_{stack}: Nx6
            
            W_A = W.A;
            W_b = W.b;
            samples = (w_stack)'; % reshape to 6 x sw_tack_len
            w_nd = size(samples,1);
            w_stack_len = size(samples,2);
            tic;

            opti = casadi.Opti( );
            beta = opti.variable( );
            y    = opti.variable(w_nd, 1);
            samplesinput = opti.parameter(w_nd, w_stack_len);
            opti.minimize(-beta);
            opti.subject_to(0 <= beta <= 1);
            opti.subject_to(W_A * y - beta * W_b <= 0);
            for i = 1:1:w_stack_len
                opti.subject_to(W_A * (samplesinput(:, i) - y) <= (1 - beta) * W_b);
            end

            opts = struct('ipopt',struct('print_level',0),'print_time',false);
            opti.solver('ipopt', opts);

            com_W_set = opti.to_function('f', {samplesinput}, {beta, y});
            tic;
            [beta_k, y_k] = com_W_set(samples);

            elapsedTime = toc;
            fprintf('coputer new W cost：%.4f s\n', elapsedTime);
            beta_k   = full(beta_k);
            y_k     = full(y_k);
            alpha_k  = 1 - beta_k;
            v_k      = y_k/beta_k;
            
            W_k = (1 - alpha_k) * Polyhedron('V', v_k') + alpha_k * W;  
end

