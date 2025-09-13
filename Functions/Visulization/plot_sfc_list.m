function plot_sfc_list(sfc_list, color)
    if nargin < 2
        color = [0.8, 0.8, 0.8];
    end
    % scatter3(sfc_list(:,1),sfc_list(:,2),sfc_list(:,3), 50, 'blue', 'filled');

    for i = 1:size(sfc_list, 1)
        row = sfc_list(i,:);
        c = row(1:3);           % center
        pos = row(4:6);         % +x, +y, +z
        neg = row(7:9);         % -x, -y, -z
        min_pt = c - neg;
        max_pt = c + pos;
        draw_box(min_pt, max_pt, color);
    end
end
