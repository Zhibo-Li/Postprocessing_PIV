%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0 = 1e-5;

cmap_1 = cmocean('ice');
cmap_2 = cmocean('solar');
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]); % plot the v-profiles

plot(nan, nan, '-', 'Color', 'k', 'LineWidth',1);
hold on
plot(nan, nan, ':', 'Color', 'k', 'LineWidth',1);
hold on

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

    U_plot_1 = U_mag(yy_start, xx_start:xx_end);
    plot(linspace(0, 1, length(U_plot_1)), U_plot_1, '-', 'Color', cmap_1((n+5)*15, :), 'LineWidth',2);
    hold on
    U_plot_2 = U_mag(yy_start:yy_end, xx_start);
    plot(linspace(0, 1, length(U_plot_2)), U_plot_2, '--', 'Color', cmap_2((n+5)*15, :), 'LineWidth',2);
    hold on
  
end

legend({'$|U|_{\rm AD}/|U_0|$', '$|U|_{\rm AB}/|U_0|$'},'FontSize', 24, ...
    'Interpreter', 'latex', 'Box','off', 'Location','northwest' )

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.3, 'FontSize', 24,'TickLabelInterpreter','latex')

ylim([0 5.0001])
ylabel('$|U|/|U_0|$','FontSize', 24,'Interpreter', 'latex');
xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\Zhibo PhD thesis\Figures\5-flexible_fiber_array' ...
    '\flow field\U_mag_profiles.pdf']);

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\Zhibo PhD thesis\' ...
    'Figures\5-flexible_fiber_array\flow field\U_mag_profiles.eps']);

% tilt_angle = 0:5:45;
% save('D:\Dropbox\GitHub\Postprocessing_PIV\Velocity_on_lattice_edge.mat', ...
%     "U_mag_up", "U_mag_left", "U_mag_max_up", "U_mag_max_left", "Pos_up", ...
%     "Pos_left", "tilt_angle")



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
ylim([-0.035 5.0035])
ylabel('$|U|_{\rm AD,\,max}/|U_0|$','FontSize', 24,'Interpreter', 'latex');

xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');
xlim([0 45.00001]); xticks([0:15:45]);

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24,'TickLabelInterpreter','latex')

ax = gca;
ax.YAxis(1).Color = [64 57 144]/255;
ax.YAxis(2).Color = [207 67 62]/255;

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\Zhibo PhD thesis\Figures\5-flexible_fiber_array' ...
    '\flow field\U_mag_max_on_lattice_edge.pdf']);



%% Plot max velocity magnitude ratio between the lattice up and left edge (to be explained)
clear; close all; clc

load('D:\Dropbox\GitHub\Postprocessing_PIV\Velocity_on_lattice_edge.mat');
for ii = 1:10

%     integral_ratio(ii) = trapz(Pos_up{ii}, U_mag_up{ii}) / trapz(Pos_left{ii}, U_mag_left{ii});
    max_Umag_ratio(ii) = U_mag_max_up(ii) / U_mag_max_left(ii);

end
% figure; plot(integral_ratio);

figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
plot(tilt_angle, max_Umag_ratio, 'o','MarkerSize', 12,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', 'm')
hold on; plot(0:45, tand(0:45), '--k','LineWidth',1)

ylim([0 1.0035])
ylabel('$|U|_{\rm AD,\,max}/|U|_{\rm AB,\,max}$','FontSize', 24,'Interpreter', 'latex');

xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');
xlim([0 45.00001]); xticks([0:15:45]);

legend({'', '$\tan\alpha_0$'},'FontSize', 24,'Interpreter', 'latex', 'Box','off', 'Location','northwest' )

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.3, 'FontSize', 24,'TickLabelInterpreter','latex')

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\Zhibo PhD thesis\Figures\5-flexible_fiber_array' ...
    '\flow field\U_mag_max_ratio.pdf']);


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

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24,'TickLabelInterpreter','latex')

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\Zhibo PhD thesis\Figures\5-flexible_fiber_array' ...
    '\flow field\Max_flowStrength_in_lattice.pdf']);



%%
% flowStrength (average)
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

    % make a mask of pillar array 
    pillar_centers = [-3e-4, 3e-4; 0, 3e-4; 3e-4, 3e-4;
        -3e-4, 0; 0, 0; 3e-4, 0;
        -3e-4, -3e-4; 0, -3e-4; 3e-4, -3e-4];
    for jj = 1:size(pillar_centers, 1)
        xx_start = find(abs(xx(1, :) - pillar_centers(jj, 1)) < 1e-8);
        yy_start = find(abs(yy(:, 1) - pillar_centers(jj, 2)) < 1e-8);
        pillar_centers_index_in_grid(jj, :) = [xx_start yy_start];
    end
    radii = 1e-4 / mean(diff(xx(1, :))); 

    mask = createCirclesMask(size(xx), pillar_centers_index_in_grid, ...
        radii*ones(size(pillar_centers_index_in_grid,1), 1));

    mask_left_edge_ind = find(abs(xx(1, :) - min(pillar_centers(:, 1))) < 1e-8);
    mask_up_edge_ind = find(abs(yy(:, 1) - min(pillar_centers(:, 2))) < 1e-8);
    mask_right_edge_ind = find(abs(xx(1, :) - max(pillar_centers(:, 1))) < 1e-8);
    mask_down_edge_ind = find(abs(yy(:, 1) - max(pillar_centers(:, 2))) < 1e-8);

    mask(:, 1:mask_left_edge_ind-1) = 0;
    mask(1:mask_up_edge_ind-1, :) = 0;
    mask(:, mask_right_edge_ind+1:end) = 0;
    mask(mask_down_edge_ind+1:end, :) = 0;

    mask = double(mask);
    mask(mask == 0) = NaN;
    flowStrength_interp_average(n) = mean(flowStrength_interp(:).*mask(:), 'omitnan')/4;

    n = n + 1;

end

tilt_angle = 0:5:45;
save('D:\Dropbox\GitHub\Postprocessing_PIV\Average_flowStrength_in_lattice.mat', ...
    "flowStrength_interp_average", "tilt_angle")



%% Plot average flowStrength in the lattice
clear; close all; clc

load('D:\Dropbox\GitHub\Postprocessing_PIV\Average_flowStrength_in_lattice.mat');

figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);

plot(tilt_angle, flowStrength_interp_average, '*', 'MarkerSize', 15, 'MarkerEdgeColor', 'k'); 
xlim([0 45]); 
ylabel('$\sigma_{\rm sum}$','FontSize', 24,'Interpreter', 'latex');
xlabel('$\alpha_0\,(^\circ)$','FontSize', 24,'Interpreter', 'latex');
xticks([0:15:45]);

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24,'TickLabelInterpreter','latex')

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\Zhibo PhD thesis\Figures\5-flexible_fiber_array' ...
    '\flow field\Average_flowStrength_in_lattice.pdf']);




function mask = createCirclesMask(varargin)
% (xDim,yDim,centers,radii)
% Create a binary mask from circle centers and radii
%
% SYNTAX:
% mask = createCirclesMask([xDim,yDim],centers,radii);
% OR
% mask = createCirclesMask(I,centers,radii);
%
% INPUTS:
% [XDIM, YDIM]   A 1x2 vector indicating the size of the desired
%                mask, as returned by [xDim,yDim,~] = size(img);
%
% I              As an alternate to specifying the size of the mask
%                (as above), you may specify an input image, I,  from which
%                size metrics are to be determined.
%
% CENTERS        An m x 2 vector of [x, y] coordinates of circle centers
%
% RADII          An m x 1 vector of circle radii
%
% OUTPUTS:
% MASK           A logical mask of size [xDim,yDim], true where the circles
%                are indicated, false elsewhere.
%
%%% EXAMPLE 1:
%   img = imread('coins.png');
%   [centers,radii] = imfindcircles(img,[20 30],...
%      'Sensitivity',0.8500,...
%      'EdgeThreshold',0.30,...
%      'Method','PhaseCode',...
%      'ObjectPolarity','Bright');
%   figure
%   subplot(1,2,1);
%   imshow(img)
%   mask = createCirclesMask(img,centers,radii);
%   subplot(1,2,2);
%   imshow(mask)
%
%%% EXAMPLE 2:
%   % Note: Mask creation is the same as in Example 1, but the image is
%   % passed in, rather than the size of the image.
%
%   img = imread('coins.png');
%   [centers,radii] = imfindcircles(img,[20 30],...
%      'Sensitivity',0.8500,...
%      'EdgeThreshold',0.30,...
%      'Method','PhaseCode',...
%      'ObjectPolarity','Bright');
%   mask = createCirclesMask(size(img),centers,radii);
%
% See Also: imfindcircles, viscircles, CircleFinder
%
% Brett Shoelson, PhD
% 9/22/2014
% Comments, suggestions welcome: brett.shoelson@mathworks.com
% Copyright 2014 The MathWorks, Inc.
narginchk(3,3)
if numel(varargin{1}) == 2
    % SIZE specified
    xDim = varargin{1}(1);
    yDim = varargin{1}(2);
else
    % IMAGE specified
    [xDim,yDim] = size(varargin{1});
end
centers = varargin{2};
radii = varargin{3};
xc = centers(:,1);
yc = centers(:,2);
[xx,yy] = meshgrid(1:yDim,1:xDim);
mask = false(xDim,yDim);
for ii = 1:numel(radii)
    mask = mask | hypot(xx - xc(ii), yy - yc(ii)) <= radii(ii);
end
end