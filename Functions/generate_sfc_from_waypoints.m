function sfc_list = generate_sfc_from_waypoints(waypoints, octomap, max_side, step, mapsize)
    % waypoints: Nx3 array of center points
    % octomap: binary occupancy map3D
    % max_side: max expansion in any direction (scalar)
    % step: expansion step size (scalar)
    % mapsize: [xmin xmax; ymin ymax; zmin zmax] – 3x2 matrix

    sfc_list = [];
    for i = 1:size(waypoints,1)
        center = waypoints(i,:);
        % Starting from step, expand the edge length evenly
        len = step;
        half_len = len / 2;
        while len <= max_side
            % generate all points within a cube, with a step size of step
            [X, Y, Z] = ndgrid(-half_len:step:half_len, -half_len:step:half_len, -half_len:step:half_len);
            cube_points = [X(:), Y(:), Z(:)] + center;      
            % Obstacle Detection
            if any(checkOccupancy(octomap, cube_points) > 0)
                len = len - step;
                half_len = len / 2;
                break; 
            end
            len = len + step;
            half_len = len / 2;
        end
        % Step 2: Based on the largest cube, expand along the y-axis
        pos_y = half_len;
        neg_y = half_len;
        while pos_y + step <= max_side/2
        
        % Construct the current cuboid point cloud, extend along the y-axis by step
        [X, Y, Z] = ndgrid(-half_len:step:half_len, -neg_y:step:pos_y+step, -half_len:step:half_len);
        cube_points = [X(:), Y(:), Z(:)] + center;
        if any(checkOccupancy(octomap, cube_points) > 0)
            pos_y = pos_y - step;
            break;
        end
        pos_y = pos_y + step;
        end

        % 继续膨胀y轴负方向
        while neg_y + step <= max_side/2
        [X, Y, Z] = ndgrid(-half_len:step:half_len, -neg_y-step:step:pos_y, -half_len:step:half_len);
        cube_points = [X(:), Y(:), Z(:)] + center;
    
        if any(checkOccupancy(octomap, cube_points) > 0)
            neg_y = neg_y - step;
            break;
        end
        neg_y = neg_y + step;
        end


        % Step 3: expand along the x-axis
        pos_x = half_len;
        neg_x = half_len;
        while pos_x + step <= max_side/2
        % Construct the current cuboid point cloud, extend along the y-axis by step
        [X, Y, Z] = ndgrid(-neg_x:step:pos_x+step, -neg_y:step:pos_y, -half_len:step:half_len);
        cube_points = [X(:), Y(:), Z(:)] + center;
        if any(checkOccupancy(octomap, cube_points) > 0)
            pos_x = pos_x - step;
            break;
        end
        pos_x = pos_x + step;
        end

        % 继续膨胀x轴负方向
        while neg_x + step <= max_side/2
        [X, Y, Z] = ndgrid(-neg_x-step:step:pos_x, -neg_y:step:pos_y, -half_len:step:half_len);
        cube_points = [X(:), Y(:), Z(:)] + center;
    
        if any(checkOccupancy(octomap, cube_points) > 0)
            neg_x = neg_x - step;
            break;
        end
        neg_x = neg_x + step;
        end

        % expand along the z-axis
        pos_z = half_len;
        neg_z = half_len;
        while pos_z + step <= max_side/2
        % Construct the current cuboid point cloud, extend along the y-axis by step
        [X, Y, Z] = ndgrid(-neg_x:step:pos_x, -neg_y:step:pos_y, -neg_z:step:pos_z+step);
        cube_points = [X(:), Y(:), Z(:)] + center;
        if any(checkOccupancy(octomap, cube_points) > 0)
            pos_z = pos_z - step;
            break;
        end
        pos_z = pos_z + step;
        end

        % 继续膨胀z轴负方向
        while neg_z + step <= max_side/2
        [X, Y, Z] = ndgrid(-neg_x:step:pos_x, -neg_y:step:pos_y, -neg_z-step:step:pos_z);
        cube_points = [X(:), Y(:), Z(:)] + center;
    
        if any(checkOccupancy(octomap, cube_points) > 0)
            neg_z = neg_z - step;
            break;
        end
        neg_z = neg_z + step;
        end
        neg_extent = [neg_x, neg_y, neg_z];
        pos_extent = [pos_x, pos_y, pos_z];
        %check the mapsize constraint
        min_pt = center - neg_extent;  % [x_min, y_min, z_min]
        max_pt = center + pos_extent;  % [x_max, y_max, z_max]
        % mapsize格式: [xmin xmax; ymin ymax; zmin zmax]
        if any(min_pt < [mapsize(1,1), mapsize(2,1), mapsize(3,1)]) || ...
           any(max_pt > [mapsize(1,2), mapsize(2,2), mapsize(3,2)])
            % 超出边界，做相应处理（如裁剪、缩小膨胀量等）
            % 例如：限制边界
            min_pt = max(min_pt, [mapsize(1,1)+0.05, mapsize(2,1)+0.05, mapsize(3,1)+0.05]);
            max_pt = min(max_pt, [mapsize(1,2)-0.05, mapsize(2,2)-0.05, mapsize(3,2)-0.05]);
            % 重新计算膨胀距离
            neg_extent = center - min_pt;
            pos_extent = max_pt - center;
        end
        
        % 最终保存
        sfc_list = [sfc_list; center, pos_extent, neg_extent];
        
    end
end
