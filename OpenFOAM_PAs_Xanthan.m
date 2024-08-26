%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0_Xanthan = 5e-3;
U_0 = 1e-5;


% U mag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Array_angle = [0 30 45]

% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
    data = readmatrix(['F:\Simulation\Xanthan\Xanthan1000ppm_Deg',num2str(Array_angle),'_PressureDiff_0o5Crop_Mid-Z0.csv']);
    Ux = data(1:end, 2);
    Uy = data(1:end, 3);
    % Uz = data(1:end, 4);
    XX = data(1:end, 5);
    YY = data(1:end, 6);
    % ZZ = data(1:end, 7);

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
    Grid_x1 = -5e-4; Grid_x2 = 5e-4; Grid_y1 = -4.5e-4; Grid_y2 = 5.5e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    Ux_interp_Xanthan = interpolant_Ux(xx,yy);
    Uy_interp_Xanthan = interpolant_Uy(xx,yy);
    

% % % % % % % % % % % % % % % % % To plot % % % % % % % % % % % % % % % %
    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    % Plot contour: |U|
    U_mag_Xanthan = sqrt(Ux_interp_Xanthan.^2+Uy_interp_Xanthan.^2);
    contourf(xx,yy,U_mag_Xanthan/max(max(U_mag_Xanthan)) ,100,'LineStyle','none');
    shading interp
    axis equal; axis off

    % Plot pillars
    hold on
%     viscircles([0 0 3 3; 0 3 0 3]'*1e-4, [1 1 1 1]*1e-4,'Color','w','LineStyle','--');
    rectangle('Position',[-1 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-1 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');

    xlim([0 3]*1e-4)
    ylim([0 3]*1e-4)

    cmocean('speed');
    caxis([0 1])
%     c = colorbar;
%     c.Label.String = '$U/U_0$';
%     c.Label.Interpreter = 'LaTeX';
%     c.TickLabelInterpreter = 'LaTeX';
%     c.FontSize = 18;

    % Plot streanlines
    hold on
    startX_ver = zeros(1, 20);
    startY_ver = linspace(1e-4,2e-4,20);
    lineobj = streamline(xx,yy,Ux_interp_Xanthan, Uy_interp_Xanthan, startX_ver(2:end-1), startY_ver(2:end-1));
    for ii = 1:numel(startX_ver)-2
        lineobj(ii).LineWidth = 0.2;
        lineobj(ii).Color = 'k';
    end

    if Array_angle ~= 0
        hold on
        startX_hor = linspace(1e-4,2e-4,20);
        startY_hor = 3e-4*ones(1, 20);
        lineobj = streamline(xx,yy,Ux_interp_Xanthan, Uy_interp_Xanthan, startX_hor(2:end-1), startY_hor(2:end-1));
        for ii = 1:numel(startX_hor)-2
            lineobj(ii).LineWidth = 0.2;
            lineobj(ii).Color = 'k';
        end
    end

    f=gcf;
    savefig(f,['F:\Simulation\Figures\Xanthan1000ppm PressureDiff' ...
        '\OpenFOAM_PAs_Xanthan_PressureDiff_Umag_Deg', num2str(Array_angle), '.fig'])
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Simulation\Figures\' ...
        'Xanthan1000ppm PressureDiff\OpenFOAM_PAs_Xanthan_PressureDiff_Umag_Deg', num2str(Array_angle), '.eps'])


    close

end


%%
% U mag difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

for Array_angle = [0 30 45]

% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
    data = readmatrix(['F:\Simulation\Xanthan\Xanthan1000ppm_Deg',num2str(Array_angle),'_PressureDiff_0o5Crop_Mid-Z0.csv']);
    Ux = data(1:end, 2);
    Uy = data(1:end, 3);
    % Uz = data(1:end, 4);
    XX = data(1:end, 5);
    YY = data(1:end, 6);
    % ZZ = data(1:end, 7);

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
    Grid_x1 = -5e-4; Grid_x2 = 5e-4; Grid_y1 = -4.5e-4; Grid_y2 = 5.5e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    Ux_interp_Xanthan = interpolant_Ux(xx,yy);
    Uy_interp_Xanthan = interpolant_Uy(xx,yy);
    

% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
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

    % Interpolate
    Ux_interp = interpolant_Ux(xx,yy);
    Uy_interp = interpolant_Uy(xx,yy);

% % % % % % % % % % % % % % % % % To plot % % % % % % % % % % % % % % % %
    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    % Plot contour: |U|
    U_mag_Xanthan = sqrt(Ux_interp_Xanthan.^2+Uy_interp_Xanthan.^2);
    U_mag = sqrt(Ux_interp.^2+Uy_interp.^2);
    contourf(xx,yy,abs(U_mag_Xanthan/max(max(U_mag_Xanthan)) ...
            - U_mag/max(max(U_mag))),100,'LineStyle','none');
    shading interp
    axis equal; axis off

    % Plot pillars
    hold on
%     viscircles([0 0 3 3; 0 3 0 3]'*1e-4, [1 1 1 1]*1e-4,'Color','w','LineStyle','--');
    rectangle('Position',[-1 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-1 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');

    xlim([0 3]*1e-4)
    ylim([0 3]*1e-4)

    cmocean('speed');
    caxis([0 1])
%     c = colorbar;
%     c.Label.String = '$U/U_0$';
%     c.Label.Interpreter = 'LaTeX';
%     c.TickLabelInterpreter = 'LaTeX';
%     c.FontSize = 18;

    f=gcf;
    savefig(f,['F:\Simulation\Figures\Xanthan1000ppm PressureDiff' ...
        '\OpenFOAM_PAs_Umag_diff_PressureDiffCase_Deg', num2str(Array_angle), '.fig'])
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Simulation\Figures\' ...
        'Xanthan1000ppm PressureDiff\OpenFOAM_PAs_Umag_diff_PressureDiffCase_Deg', num2str(Array_angle), '.eps'])

    close

end



%%
% flowType
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

for Array_angle = [0 30 45]

% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
    data = readmatrix(['F:\Simulation\Xanthan\Xanthan1000ppm_Deg',num2str(Array_angle),'_PressureDiff_0o5Crop_Mid-Z0.csv']);
    flowType = data(1:end, 1);
    XX = data(1:end, 5);
    YY = data(1:end, 6);
    % ZZ = data(1:end, 7);

    % Rotation matrix
    RotMatrix = rotz(-Array_angle); RotMatrix = RotMatrix(1:2, 1:2);

    % Rotate the position
    Rotated_pos = RotMatrix * [XX,YY]';
    Rotated_xx = Rotated_pos(1, :);
    Rotated_yy = Rotated_pos(2, :);

    % Interpolate to grid
    interpolant_flowType = scatteredInterpolant(Rotated_xx',Rotated_yy',flowType,'natural','none');

    % Grid
    Grid_density = 1001;
    Grid_x1 = -5e-4; Grid_x2 = 5e-4; Grid_y1 = -4.5e-4; Grid_y2 = 5.5e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    flowType_interp_Xanthan = interpolant_flowType(xx,yy);


% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
    data = readmatrix(['Z:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\Data\Deg', num2str(Array_angle), '_0o5Crop_Mid-Z0.csv']);
    flowType = data(1:end, 2);
    XX = data(1:end, 6);
    YY = data(1:end, 7);
    % ZZ = data(1:end, 8);

    % Rotation matrix
    RotMatrix = rotz(-Array_angle); RotMatrix = RotMatrix(1:2, 1:2);

    % Rotate the position
    Rotated_pos = RotMatrix * [XX,YY]';
    Rotated_xx = Rotated_pos(1, :);
    Rotated_yy = Rotated_pos(2, :);

    % Interpolate to grid
    interpolant_flowType = scatteredInterpolant(Rotated_xx',Rotated_yy',flowType,'natural','none');

    % Interpolate
    flowType_interp = interpolant_flowType(xx,yy);
        
% % % % % % % % % % % % % % % % % To plot % % % % % % % % % % % % % % % %
    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    % Plot contour: flowType
    contourf(xx,yy,flowType_interp_Xanthan,100,'LineStyle','none');
    shading interp
    axis equal; axis off

    % Plot pillars
    hold on
%     viscircles([0 0 3 3; 0 3 0 3]'*1e-4, [1 1 1 1]*1e-4,'Color','w','LineStyle','--');
    rectangle('Position',[-1 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-1 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');

    xlim([0 3]*1e-4)
    ylim([0 3]*1e-4)

    cmocean('balance');
    caxis([-1 1])
    % c = colorbar;
    % c.Label.String = '$flowType$';
    % c.Label.Interpreter = 'LaTeX';
    % c.TickLabelInterpreter = 'LaTeX';
    % c.FontSize = 18;

    f=gcf;
    savefig(f,['F:\Simulation\Figures\Xanthan1000ppm PressureDiff' ...
        '\OpenFOAM_PAs_Xanthan_PressureDiff_flowType_Deg', num2str(Array_angle), '.fig'])
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Simulation\Figures\Xanthan1000ppm PressureDiff' ...
        '\OpenFOAM_PAs_Xanthan_PressureDiff_flowType_Deg', num2str(Array_angle), '.eps'])

    close

end