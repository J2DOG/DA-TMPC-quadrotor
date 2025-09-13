function flag = box_is_overlap(a, b)
    % Determine if b intersects with a
    ca = a(1:3); pa = a(4:6); na = a(7:9);
    cb = b(1:3); pb = b(4:6); nb = b(7:9);

    min_a = ca - na;
    max_a = ca + pa;
    min_b = cb - nb;
    max_b = cb + pb;

    % Check that all three-dimensional intervals have intersections
    flag = all(min_a <= max_b) && all(max_a >= min_b);
end