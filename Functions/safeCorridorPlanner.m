classdef safeCorridorPlanner < handle
    % SAFECORRIDORPLANNER Summary of this class goes here
    % Detailed explanation goes here
    
    properties
        obstacle_list
        obstacle_num
        map_poly
        ini_global_path
        d_min
        expandsolver
    end
    
    methods
        function obj = safeCorridorPlanner(obstacle_list, map_poly, ini_global_path, d_min)
            %SAFECORRIDORPLANNER Construct an instance of this class 
            obj.obstacle_list = obstacle_list;
            obj.obstacle_num = length(obstacle_list);
            obj.map_poly = map_poly;
            obj.ini_global_path = ini_global_path;
            obj.d_min = d_min;
            obj.expandsolver = obj.build_expansion_solver(size(map_poly.A,1), size(map_poly.A,2), obstacle_list, d_min);
        end
        
        function corridors = plan(obj)
            ini_path_length = size(obj.ini_global_path,1);
            corridors = Polyhedron.empty(); % polyhedron array
            k = 0;
            i = 1;
            max_inner_iter = 100; 
        
            while i < ini_path_length
                k = k + 1;
                p_curr = obj.ini_global_path(i,:)'; 
        
                alpha_opt = obj.getAlphaOpt(p_curr);
        
                % --- Construct the Half-spce representation of Safe set: A x <= alpha*b + (1-alpha)*A*p_curr
                A_map = obj.map_poly.A;
                b_map = obj.map_poly.b;
                SFC_b = alpha_opt * b_map + (1 - alpha_opt) * (A_map * p_curr);
                SFC = Polyhedron('A', A_map, 'b', SFC_b);
        
                fprintf('SFC No.%d, i = %d,  alpha_opt = %.4f, point = [%.2f %.2f %.2f]\n', ...
                        k,i, alpha_opt, p_curr(1), p_curr(2), p_curr(3));
                corridors(end+1) = SFC;
        
                % If at the end point, exit
                if i >= ini_path_length
                    break;
                end
        
                % Check if next waypoint is inner the safe set （Ax <= b）
                p_next = obj.ini_global_path(i+1,:)';
                violation = SFC.A * p_next - SFC.b;
        
                % If the next point is outside, try to find the insert point
                % insert point is the line segment from waypoint now to next waypoint intersect with the safe set
                % and construct a new Safe set at the intersection point.
                inner_iter = 0;
                while any(violation > 1e-9)  
                    inner_iter = inner_iter + 1;
                    if inner_iter > max_inner_iter
                        warning('Inner loop reached max iterations; breaking to avoid infinite loop.');
                        break;
                    end
        
                    L = Polyhedron('V', [p_curr'; p_next']);  % 1D polyhedron (segment)
                    Pint = intersect(SFC, L);
                    V = Pint.V; % n x 3
                    if isempty(V)
                        warning('Intersection polyhedron has no vertices (unexpected).');
                        break;
                    end
                    dists = vecnorm((V.' - p_curr), 2, 1); % 1 x n
                    [~, idx] = max(dists);

                    p_insert = V(idx, :)';
                    p_curr = p_insert;
                    k = k + 1;
                    alpha_opt = obj.getAlphaOpt(p_curr);
                    SFC_b = alpha_opt * b_map + (1 - alpha_opt) * (A_map * p_curr);
                    SFC = Polyhedron('A', A_map, 'b', SFC_b);
    
                    fprintf('InsertPoint:\n');
                    fprintf('SFC No.%d, i = %d, alpha_opt = %.4f, point = [%.2f %.2f %.2f], point_next = [%.2f %.2f %.2f]\n', ...
                            k,i, alpha_opt, p_curr(1), p_curr(2), p_curr(3), p_next(1), p_next(2), p_next(3));
                    corridors(end+1) = SFC;
    
                    violation = SFC.A * p_next - SFC.b;
                    

                end 
        
                while i < ini_path_length
                    p_next = obj.ini_global_path(i+1,:)';
                    violation = SFC.A * p_next - SFC.b;
                    if any(violation > 1e-9)
                        break; 
                    end
                    i = i + 1; 
                end
            end 
        end

        function [corr_simple, sequence] = planSimpled(obj)
            % Step 1: simplify corridors
            corridors = obj.plan();   % 1xN mpt3 Polyhedron
        
            while true
                removed = false;
                for i = 2:length(corridors)-1
                    intersec = intersect(corridors(i-1), corridors(i+1));
                    if ~intersec.isEmptySet
                        corridors(i) = [];
                        removed = true;
                        break;  % Delete one and start over from scratch.
                    end
                end
                if ~removed
                    break;
                end
            end
        
            corr_simple = corridors;
        
            sequence = [];
            for i = 1:length(corridors)-1
                intersec = intersect(corridors(i), corridors(i+1));
                if ~intersec.isEmptySet
                    center = mean(intersec.V,1);  
                    sequence = [sequence; center];
                end
            end
            sequence = [obj.ini_global_path(1,:); sequence; obj.ini_global_path(end,:)];
        end

        function [corr_expand] = planExpand(obj)
            corridors = obj.planSimpled();
            for i = 1:length(corridors)
                corr_expand = obj.expand(corridors);
            end
        end
        function corridors_expand = expand(obj, corridors)
            % Still have bug, don't use this method!!!

            
            n_corridors = length(corridors);
            corridors_expand = repmat(Polyhedron(), size(corridors));
            
            param_base = {};
            for i = 1:obj.obstacle_num
                A_obs = obj.obstacle_list{i}.A;
                b_obs = obj.obstacle_list{i}.b;
                param_base{end+1} = A_obs;
                param_base{end+1} = b_obs;
            end
            

            for i = 1:n_corridors
                P = corridors(i);
                S_val = P.A;  % Sx <= s
                s_val = P.b;
                input_val = [{S_val, s_val}, param_base{:}];
            
                alpha_opt = obj.expandsolver(input_val{:});
                s_new = (1 + alpha_opt{1}) .* s_val;
                s_new = full(s_new);
                corridors_expand(i) = Polyhedron('A', S_val, 'b', s_new);
            end
    end

    function expandsolver = build_expansion_solver(obj,m, n, obstacle_list, d_min)
            % still unfinished
            
            import casadi.*
            
            opti = casadi.Opti();
            
            % 参数：Sx <= s
            S = opti.parameter(m, n);
            s = opti.parameter(m, 1);
            M_map = obj.map_poly.A;
            m_map = obj.map_poly.b;
            params = {S, s};
            

            alpha = opti.variable(m, 1);
            
            opti.minimize(-sum(alpha));
            [p, ~] = size(M_map);    
            H = opti.variable(p, m);
            s_expand = (1 + alpha) .* s;
            
            opti.subject_to(H * S == M_map);              % 保证法向组合正确
            opti.subject_to(H * s_expand <= m_map);   % 偏置约束不能超出 map
            opti.subject_to(H(:) >= 0);                  % H 元素非负
            for i = 1:length(obstacle_list)
                A_obs = obstacle_list{i}.A;
                b_obs = obstacle_list{i}.b;
                li = size(A_obs, 1);
            
                Oi = opti.parameter(li, n);
                oi = opti.parameter(li, 1);
                params = [params, {Oi, oi}];
            
                mu = opti.variable(m, 1);
                lambda = opti.variable(li, 1);
            
                support_term = (1 + alpha) .* s;
                dis = -support_term' * mu - oi' * lambda;
                opti.subject_to(dis >= d_min);
                opti.subject_to(S' * mu + Oi' * lambda == 0);
                opti.subject_to(sumsqr(Oi' * lambda) <= 1);
                opti.subject_to(mu >= 0);
                opti.subject_to(lambda >= 0);
            end
            
            opti.subject_to(alpha >= 0);  
            opts = struct;
            opts.ipopt.print_level = 0;        
            opts.print_time = false;           
            opti.solver('ipopt', opts);
            
            expandsolver = opti.to_function('expansion_solver', params, {alpha});
            end



        function alpha_opt = getAlphaOpt(obj,y_s)
            alpha_up = 1;
            alpha_down = 0;
            error_tolerance = 1e-4;
            while alpha_up - alpha_down > error_tolerance
                alpha_mid = 0.5 * (alpha_down + alpha_up);
                if obj.checkAlphaFeasible(y_s, alpha_mid)
                    alpha_down = alpha_mid;
                else
                    alpha_up = alpha_mid;
                end
            end
            alpha_opt = alpha_down;
        end
             

        function flag = checkAlphaFeasible(obj, y_s_val, alpha_val)
            import casadi.*
            Fy = obj.map_poly.A;
            fy = obj.map_poly.b;
            opti  = casadi.Opti();
        
            y_s = opti.parameter(size(Fy,2),1);
            alpha = opti.parameter(1,1);

            lambda = cell(obj.obstacle_num, 1);
            mui = cell(obj.obstacle_num, 1);
    
            opti.minimize(0);

            for i  = 1:obj.obstacle_num
                mui{i} = opti.variable(size(fy, 1), 1);
                lambda{i} = opti.variable(size(obj.obstacle_list{i}.A, 1), 1);
                dis = -(alpha*fy)'*mui{i} + (obj.obstacle_list{i}.A*(1-alpha)*y_s-obj.obstacle_list{i}.b)'*lambda{i};
                opti.subject_to(dis >= obj.d_min);
                opti.subject_to(Fy'*mui{i}+obj.obstacle_list{i}.A'*lambda{i} == 0);
                opti.subject_to(sumsqr(obj.obstacle_list{i}.A'*lambda{i}) <= 1);
                opti.subject_to(mui{i} >= 0);
                opti.subject_to(lambda{i} >= 0);
            end
            opts = struct('ipopt', struct('print_level', 0), 'print_time', false);
            opti.solver('ipopt', opts);
        
            opti.set_value(y_s, y_s_val);
            opti.set_value(alpha, alpha_val);
            try
                sol = opti.solve();        
                flag = 1;              
            catch
                flag = 0; 
            end 
        end

    end
end

