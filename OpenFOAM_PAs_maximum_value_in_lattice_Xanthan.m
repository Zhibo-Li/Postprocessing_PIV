%% Xanthan (flow rate driven) and Newtonian fluid
%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0_Xanthan = 5e-3;
U_0 = 1e-5;

cmap_1 = cmocean('ice');
cmap_2 = cmocean('solar');
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]); % plot the v-profiles

% % % plot(nan, nan, '-', 'Color', 'k', 'LineWidth',1);
% % % hold on
% % % plot(nan, nan, ':', 'Color', 'k', 'LineWidth',1);
% % % hold on

n = 1;
for Array_angle = [0 30 45]

% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
    data = readmatrix(['F:\Simulation\Xanthan\Xanthan1000ppm_Deg',num2str(Array_angle),'_0o5Crop_Mid-Z0.csv']);
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

    % U_mag
    U_mag_Xanthan = sqrt(Ux_interp_Xanthan.^2+Uy_interp_Xanthan.^2)/U_0_Xanthan;
    

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

    % U_mag
    U_mag = sqrt(Ux_interp.^2+Uy_interp.^2)/U_0;

    xx_start = find(abs(xx(1, :) - 0) < 1e-8); 
    xx_end = find(abs(xx(1, :) - 3e-4) < 1e-8);
    yy_start = find(abs(yy(:, 1) - 0) < 1e-8); 
    yy_end = find(abs(yy(:, 1) - 3e-4) < 1e-8);

    U_plot_1 = U_mag(yy_start, xx_start:xx_end);
    plot(linspace(0, 1, length(U_plot_1)), U_plot_1, '-', 'Color', cmap_1(n*70, :), 'LineWidth',2);
    hold on
    U_plot_3 = U_mag_Xanthan(yy_start, xx_start:xx_end);
    plot(linspace(0, 1, length(U_plot_3)), U_plot_3, ':', 'Color', cmap_1(n*70, :), 'LineWidth',2);
    hold on
    U_plot_2 = U_mag(yy_start:yy_end, xx_start);
    plot(linspace(0, 1, length(U_plot_2)), U_plot_2, '-', 'Color', cmap_2(n*70, :), 'LineWidth',2);
    hold on
    U_plot_4 = U_mag_Xanthan(yy_start:yy_end, xx_start);
    plot(linspace(0, 1, length(U_plot_4)), U_plot_4, ':', 'Color', cmap_2(n*70, :), 'LineWidth',2);
    hold on

    n = n + 1;
  
end

legend({'AD ($0^\circ$, Newtonian)', 'AD ($0^\circ$, Xanthan)', ...
    'AB ($0^\circ$, Newtonian)', 'AB ($0^\circ$, Xanthan)', ...
    'AD ($30^\circ$, Newtonian)', 'AD ($30^\circ$, Xanthan)', ...
    'AB ($30^\circ$, Newtonian)', 'AB ($30^\circ$, Xanthan)', ...
    'AD ($45^\circ$, Newtonian)', 'AD ($45^\circ$, Xanthan)', ...
    'AB ($45^\circ$, Newtonian)', 'AB ($45^\circ$, Xanthan)'},'FontSize', 18, ...
    'Interpreter', 'latex', 'Box','off', 'Location','bestoutside' )

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.3, 'FontSize', 24,'TickLabelInterpreter','latex')

ylim([0 5.0001])
ylabel('$|U|/|U_0|$','FontSize', 24,'Interpreter', 'latex');
% xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Simulation\Figures\U_mag_profiles.pdf']);

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Simulation\Figures\U_mag_profiles.eps']);



%% Xanthan: pressure driven and flow rate driven
%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0_Xanthan = 5e-3;
U_0 = 1e-5;

cmap_1 = cmocean('ice');
cmap_2 = cmocean('solar');
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]); % plot the v-profiles

% % % plot(nan, nan, '-', 'Color', 'k', 'LineWidth',1);
% % % hold on
% % % plot(nan, nan, ':', 'Color', 'k', 'LineWidth',1);
% % % hold on

n = 1;
for Array_angle = [0 30 45]

% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
    data = readmatrix(['F:\Simulation\Xanthan\Xanthan1000ppm_Deg',num2str(Array_angle),'_0o5Crop_Mid-Z0.csv']);
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

    % U_mag
    U_mag_Xanthan = sqrt(Ux_interp_Xanthan.^2+Uy_interp_Xanthan.^2);
    U_mag_Xanthan = U_mag_Xanthan/max(max(U_mag_Xanthan));
    

% % % % % % % % % % % % % % read simulation data % % % % % % % % % % % % %
%     data = readmatrix(['F:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
%         num2str(Array_angle), 'deg\Data\Deg', num2str(Array_angle), '_0o5Crop_Mid-Z0.csv']);
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

    % Interpolate
    Ux_interp = interpolant_Ux(xx,yy);
    Uy_interp = interpolant_Uy(xx,yy);

    % U_mag
    U_mag = sqrt(Ux_interp.^2+Uy_interp.^2);
    U_mag = U_mag/max(max(U_mag));

    xx_start = find(abs(xx(1, :) - 0) < 1e-8); 
    xx_end = find(abs(xx(1, :) - 3e-4) < 1e-8);
    yy_start = find(abs(yy(:, 1) - 0) < 1e-8); 
    yy_end = find(abs(yy(:, 1) - 3e-4) < 1e-8);

    U_plot_1 = U_mag(yy_start, xx_start:xx_end);
    plot(linspace(0, 1, length(U_plot_1)), U_plot_1, '-', 'Color', cmap_1(n*70, :), 'LineWidth',2);
    hold on
    U_plot_3 = U_mag_Xanthan(yy_start, xx_start:xx_end);
    plot(linspace(0, 1, length(U_plot_3)), U_plot_3, ':', 'Color', cmap_1(n*70, :), 'LineWidth',2);
    hold on
    U_plot_2 = U_mag(yy_start:yy_end, xx_start);
    plot(linspace(0, 1, length(U_plot_2)), U_plot_2, '-', 'Color', cmap_2(n*70, :), 'LineWidth',2);
    hold on
    U_plot_4 = U_mag_Xanthan(yy_start:yy_end, xx_start);
    plot(linspace(0, 1, length(U_plot_4)), U_plot_4, ':', 'Color', cmap_2(n*70, :), 'LineWidth',2);
    hold on

    n = n + 1;
  
end

legend({'AD ($0^\circ$, Pressure)', 'AD ($0^\circ$, Flow Rate)', ...
    'AB ($0^\circ$, Pressure)', 'AB ($0^\circ$, Flow Rate)', ...
    'AD ($30^\circ$, Pressure)', 'AD ($30^\circ$, Flow Rate)', ...
    'AB ($30^\circ$, Pressure)', 'AB ($30^\circ$, Flow Rate)', ...
    'AD ($45^\circ$, Pressure)', 'AD ($45^\circ$, Flow Rate)', ...
    'AB ($45^\circ$, Pressure)', 'AB ($45^\circ$, Flow Rate)'},'FontSize', 18, ...
    'Interpreter', 'latex', 'Box','off', 'Location','bestoutside' )

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.3, 'FontSize', 24,'TickLabelInterpreter','latex')

ylim([0 1])
ylabel('$|U|/|U_{\rm{max}}|$','FontSize', 24,'Interpreter', 'latex');
% xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Simulation\Figures\U_mag_profiles_Xanthan.pdf']);

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Simulation\Figures\U_mag_profiles_Xanthan.eps']);

