function d = optimize_translation(A,b,points)
    % By translating each hyperplane to use the input polyhedron(Ax < b) to fit the given convex set of points.
    % Transforming a convex polyhedron represented by Ax < b into Ax < b + d to bounding all the given points.
    fprintf('\n');
    fprintf('**************\n');
    fprintf(' Fitting the Set.\n');
    
    nc = size(A,1);
    num_points = size(points,1);
    % Calculate the geometric center of a point set
    centroid = mean(points, 1)';
    % Translate polyhedron from origin to the geometric center of the point set
    b = b + A*centroid;
    yalmip('clear')
    d = sdpvar(nc,1);
    constraints = [];
    for i = 1:nc
        for j = 1:num_points
            constraints = [constraints, A(i,:)*points(j,:)' <= b(i) + d(i)];
        end
    end
    % Objective function: L2 norm of d
    objective = sum(d.^2);  
    % Solve - Use a quadratic programming solver
    options = sdpsettings('solver', 'quadprog', 'verbose', 0);
    options.quadprog.MaxIterations = 10000;
    options.quadprog.OptimalityTolerance = 1e-12;
    options.quadprog.StepTolerance = 1e-12;
    
    solution = optimize(constraints, objective, options);
    if solution.problem == 0
        d = value(d);
        fprintf('Optimization success:||d||_2 = %.6f\n', norm(d, 2));
    else
        error('fail: %s', solution.info);
    end
    fprintf('[%s]  Finished.\n', datestr(now, 'HH:MM:SS'));
    fprintf('**************\n');
end