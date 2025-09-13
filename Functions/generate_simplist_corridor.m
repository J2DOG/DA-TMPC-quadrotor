function [sfc_list_simplist, overlap_box] = generate_simplist_corridor(sfc_list)
     while true
        N = size(sfc_list, 1);
        removed = false;

        for i = 2:N-1
            prev = sfc_list(i-1, :);
            curr = sfc_list(i, :);
            next = sfc_list(i+1, :);

            if box_is_overlap(prev, next)
                sfc_list(i, :) = [];  % 删除当前冗余 corridor
                removed = true;
                break;  % 删除一个就立即退出，重新从头开始检查
            end
        end

        if ~removed
            break;  % 没有再删除任何 corridor，结束循环
        end
     end
     sfc_list_simplist = sfc_list;
     
     % === Step 2: 计算每一对相邻 corridor 的重叠区域 ===
    overlap_box = [];
    for i = 1:size(sfc_list_simplist,1)-1
        c1 = sfc_list_simplist(i, 1:3);
        p1 = sfc_list_simplist(i, 4:6);
        n1 = sfc_list_simplist(i, 7:9);
        min1 = c1 - n1;
        max1 = c1 + p1;

        c2 = sfc_list_simplist(i+1, 1:3);
        p2 = sfc_list_simplist(i+1, 4:6);
        n2 = sfc_list_simplist(i+1, 7:9);
        min2 = c2 - n2;
        max2 = c2 + p2;

        % 重叠计算
        overlap_min = max(min1, min2);
        overlap_max = min(max1, max2);

        if all(overlap_min <= overlap_max)
            center = (overlap_min + overlap_max) / 2;
            pos_extent = max(overlap_max - center, 0);
            neg_extent = max(center - overlap_min, 0);
            box_overlap = [center, pos_extent, neg_extent];
            overlap_box = [overlap_box; box_overlap];
        end
    end

end

