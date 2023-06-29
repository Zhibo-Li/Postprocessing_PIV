%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0 = 1e-5;

for Array_angle = 0:5:45

    % read simulation data
    data = readmatrix(['F:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\Data\Deg', num2str(Array_angle), '_0o5Crop_Mid-Z0.csv']);
    Ux = data(1:end, 3);
    Uy = data(1:end, 4);
    % Uz = data(1:end, 5);
    XX = data(1:end, 6);
    YY = data(1:end, 7);
    % ZZ = data(1:end, 8);

    % Rotation matrix
    RotMatrix = rotz(-Array_angle); RotMatrix = RotMatrix(1:2, 1:2);

    % Rotate the velocity
    Rotated_velo = RotMatrix * [Ux, Uy]';
    Rotated_Ux = Rotated_velo(1, :);
    Rotated_Uy = Rotated_velo(2, :);

    % Rotate the position
    Rotated_pos = RotMatrix * [XX,YY]';
    Rotated_xx = Rotated_pos(1, :);
    Rotated_yy = Rotated_pos(2, :);

    % Interpolate to grid
    interpolant_Ux = scatteredInterpolant(Rotated_xx',Rotated_yy',Rotated_Ux','natural','none');
    interpolant_Uy = scatteredInterpolant(Rotated_xx',Rotated_yy',Rotated_Uy','natural','none');

    % Grid
    Grid_density = 1001;
    Grid_x1 = -6e-4; Grid_x2 = 6e-4; Grid_y1 = 0; Grid_y2 = 12e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    Ux_interp = interpolant_Ux(xx,yy);
    Uy_interp = interpolant_Uy(xx,yy);

    % To plot
    figure('color', 'w'); set(gcf, 'Position', [100 100 800 400]);
    axis equal; axis off

    % Plot pillars
    hold on
    rectangle('Position',[-1 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-1 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-4 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-4 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
%     rectangle('Position',[-4 -4 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
%     rectangle('Position',[-1 -4 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
%     rectangle('Position',[2 -4 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');

    xlim([-3 3]*1e-4)
    ylim([0 3]*1e-4)
    % Plot streanlines
    hold on
    startX_ver_1 = -4.25e-4*ones(1, 100);
    startY_ver_1 = linspace(0,12e-4,100);
    lineobj = streamline(xx,yy,Ux_interp, Uy_interp, startX_ver_1(2:end-1), startY_ver_1(2:end-1));
    for ii = 1:numel(startX_ver_1)-2
        lineobj(ii).LineWidth = 1;
        lineobj(ii).Color = 'k';
        lineobj(ii).LineStyle = ':';
    end
    hold on

    f=gcf;
    savefig(f,['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
        '\Dual lattice flow field\OpenFOAM_PAs_U_Deg', num2str(Array_angle), '_dual_lattice.fig'])
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
        'Figures\5-flexible_fiber_array\Dual lattice flow field\OpenFOAM_PAs_U_Deg', num2str(Array_angle), '_dual_lattice.eps'])

    close

end



%% threefold-lattice

clear; close all; clc

% Initial velocity
U_0 = 1e-5;

for Array_angle = 20

    % read simulation data
    data = readmatrix(['F:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\Data\Deg', num2str(Array_angle), '_0o5Crop_Mid-Z0.csv']);
    Ux = data(1:end, 3);
    Uy = data(1:end, 4);
    % Uz = data(1:end, 5);
    XX = data(1:end, 6);
    YY = data(1:end, 7);
    % ZZ = data(1:end, 8);

    % Rotation matrix
    RotMatrix = rotz(-Array_angle); RotMatrix = RotMatrix(1:2, 1:2);

    % Rotate the velocity
    Rotated_velo = RotMatrix * [Ux, Uy]';
    Rotated_Ux = Rotated_velo(1, :);
    Rotated_Uy = Rotated_velo(2, :);

    % Rotate the position
    Rotated_pos = RotMatrix * [XX,YY]';
    Rotated_xx = Rotated_pos(1, :);
    Rotated_yy = Rotated_pos(2, :);

    % Interpolate to grid
    interpolant_Ux = scatteredInterpolant(Rotated_xx',Rotated_yy',Rotated_Ux','natural','none');
    interpolant_Uy = scatteredInterpolant(Rotated_xx',Rotated_yy',Rotated_Uy','natural','none');

    % Grid
    Grid_density = 1001;
    Grid_x1 = -6e-4; Grid_x2 = 6e-4; Grid_y1 = 0; Grid_y2 = 12e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    Ux_interp = interpolant_Ux(xx,yy);
    Uy_interp = interpolant_Uy(xx,yy);

    % To plot
    figure('color', 'w'); set(gcf, 'Position', [100 100 800 400]);
    axis equal; axis off

    % Plot pillars
    hold on
    rectangle('Position',[-1 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-1 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-4 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-4 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[5 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[5 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
%     rectangle('Position',[2 -4 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');

    xlim([-3 6]*1e-4)
    ylim([0 3]*1e-4)
    % Plot streanlines
    hold on
    startX_ver_1 = -4.1e-4*ones(1, 200);
    startY_ver_1 = linspace(0,12e-4,200);
    lineobj = streamline(xx,yy,Ux_interp, Uy_interp, startX_ver_1(2:end-1), startY_ver_1(2:end-1));
    for ii = 1:numel(startX_ver_1)-2
        lineobj(ii).LineWidth = 1;
        lineobj(ii).Color = 'k';
        lineobj(ii).LineStyle = ':';
    end
    hold on

    f=gcf;
    savefig(f,['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
        '\Dual lattice flow field\OpenFOAM_PAs_U_Deg', num2str(Array_angle), '_threefold_lattice.fig'])
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
        'Figures\5-flexible_fiber_array\Dual lattice flow field\OpenFOAM_PAs_U_Deg', num2str(Array_angle), '_threefold_lattice.eps'])

    close

end
