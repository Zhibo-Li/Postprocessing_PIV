%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0 = 1e-5;

n = 1;
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
    Grid_x1 = -5e-4; Grid_x2 = 5e-4; Grid_y1 = -4.5e-4; Grid_y2 = 5.5e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    Ux_interp = interpolant_Ux(xx,yy);
    Uy_interp = interpolant_Uy(xx,yy);

    % U_mag
    U_mag = sqrt(Ux_interp.^2+Uy_interp.^2)/U_0;

    xx_start = find(abs(xx(1, :) - 0) < 1e-8); 
    xx_end = find(abs(xx(1, :) - 3e-4) < 1e-8);
    yy_start = find(abs(yy(:, 1) - 0) < 1e-8); 
    yy_end = find(abs(yy(:, 1) - 3e-4) < 1e-8);

    U_mag_max_left(n) = max(U_mag(yy_start:yy_end, xx_start));
    U_mag_max_up(n) = max(U_mag(yy_start, xx_start:xx_end));

    U_mag_left{n} = U_mag(yy_start:yy_end, xx_start);
    U_mag_up{n} = U_mag(yy_start, xx_start:xx_end);
    Pos_left{n} = xx(yy_start, xx_start:xx_end);
    Pos_up{n} = yy(yy_start:yy_end, xx_start);
    n = n + 1;

%     figure; 
%     plot(xx(yy_start, xx_start:xx_end), U_mag(yy_start, xx_start:xx_end), 'r');
%     hold on
%     plot(yy(yy_start:yy_end, xx_start), U_mag(yy_start:yy_end, xx_start), 'g');
  
end

tilt_angle = 0:5:45;
save('D:\Dropbox\GitHub\Postprocessing_PIV\Velocity_on_lattice_edge.mat', ...
    "U_mag_up", "U_mag_left", "U_mag_max_up", "U_mag_max_left", "Pos_up", ...
    "Pos_left", "tilt_angle")



%% Plot max velocity magnitude on the lattice edge (left and up)
clear; close all; clc

load('D:\Dropbox\GitHub\Postprocessing_PIV\Velocity_on_lattice_edge.mat');

figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);

yyaxis left
plot(tilt_angle, U_mag_max_left, 'o','MarkerSize', 12,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', [64 57 144]/255); 
ylim([0 5.0035])
ylabel('$|U|_{\rm AB,\,max}/|U_0|$','FontSize', 24,'Interpreter', 'latex');

yyaxis right
plot(tilt_angle, U_mag_max_up, 'o','MarkerSize', 12,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', [207 67 62]/255);
ylim([0 5.0035])
ylabel('$|U|_{\rm AD,\,max}/|U_0|$','FontSize', 24,'Interpreter', 'latex');

xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');
xlim([0 45.00001]); xticks([0:15:45]);

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

ax = gca;
ax.YAxis(1).Color = [64 57 144]/255;
ax.YAxis(2).Color = [207 67 62]/255;

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
    '\flow field\U_mag_max_on_lattice_edge.pdf']);



%% Plot max velocity magnitude ratio between the lattice up and left edge (to be explained)
clear; close all; clc

load('D:\Dropbox\GitHub\Postprocessing_PIV\Velocity_on_lattice_edge.mat');
for ii = 1:10

    integral_ratio(ii) = trapz(Pos_up{ii}, U_mag_up{ii}) / trapz(Pos_left{ii}, U_mag_left{ii});
    max_Umag_ratio(ii) = U_mag_max_up(ii) / U_mag_max_left(ii);

end
figure; plot(integral_ratio);

figure; plot(max_Umag_ratio)
hold on; plot(tand(tilt_angle))



%%
% flowStrength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

n = 1;
for Array_angle = 0:5:45

    % read simulation data
    data = readmatrix(['Z:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\Data\Deg', num2str(Array_angle), '_0o5Crop_Mid-Z0.csv']);
    flowStrength = data(1:end, 1);
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
    interpolant_flowStrength = scatteredInterpolant(Rotated_xx',Rotated_yy',flowStrength,'natural','none');

    % Grid
    Grid_density = 1001;
    Grid_x1 = -5e-4; Grid_x2 = 5e-4; Grid_y1 = -4.5e-4; Grid_y2 = 5.5e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    flowStrength_interp = interpolant_flowStrength(xx,yy);

    flowStrength_interp_max(n) = max(max(flowStrength_interp));
    n = n + 1;

end

tilt_angle = 0:5:45;
save('D:\Dropbox\GitHub\Postprocessing_PIV\Max_flowStrength_in_lattice.mat', ...
    "flowStrength_interp_max", "tilt_angle")



%% Plot max flowStrength in the lattice
clear; close all; clc

load('D:\Dropbox\GitHub\Postprocessing_PIV\Max_flowStrength_in_lattice.mat');

figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);

plot(tilt_angle, flowStrength_interp_max, '*','MarkerSize', 15,'MarkerEdgeColor','k'); 
xlim([0 45]); ylim([2 4])
ylabel('$\sigma_{\rm max}$','FontSize', 24,'Interpreter', 'latex');
xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');
xticks([0:15:45]);

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
    '\flow field\Max_flowStrength_in_lattice.pdf']);

