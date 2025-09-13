function draw_box(min_pt, max_pt, color)
    % Draw a 3D box from min_pt to max_pt

    x = [min_pt(1) max_pt(1)];
    y = [min_pt(2) max_pt(2)];
    z = [min_pt(3) max_pt(3)];

    % Corners of the box
    verts = [
        x(1) y(1) z(1);
        x(2) y(1) z(1);
        x(2) y(2) z(1);
        x(1) y(2) z(1);
        x(1) y(1) z(2);
        x(2) y(1) z(2);
        x(2) y(2) z(2);
        x(1) y(2) z(2)];

    faces = [ 1 2 3 4;   % bottom
              5 6 7 8;   % top
              1 2 6 5;   % side
              2 3 7 6;
              3 4 8 7;
              4 1 5 8];

    patch('Vertices', verts, 'Faces', faces, ...
          'FaceColor', color, 'FaceAlpha', 0.2, ...
          'EdgeColor', color, 'LineWidth', 1,'HandleVisibility','off');
end
