%%% OpenFOAM data postprocessing: microchannel with an individual obstacle
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% include the corresponding uPIV data (20221031)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% read simulation data
data = readmatrix('F:\Simulation\202302_triangular_pillar_EXPconfig\Data\MidPlane_Velocity\MidPlane_Velocity0.csv');
Ux = data(1:end, 2); 
Uy = data(1:end, 3); 
% Uz = data(1:end, 4); 
XX = data(1:end, 5);
YY = data(1:end, 6);
% ZZ = data(1:end, 7);

% read obstacle position data
obs_data = readmatrix('F:\Simulation\202302_triangular_pillar_EXPconfig\Data\MidPlane_Velocity\MidPlane_Velocity5.csv');
obs_XX = obs_data(1:end, 5);
obs_YY = obs_data(1:end, 6);
% sort the coordinates clockwise
obs_center_XX = mean(obs_XX); obs_center_YY = mean(obs_YY); 
[theta, ~] = cart2pol(obs_XX-obs_center_XX, obs_YY-obs_center_YY);
obs_2d = [obs_XX, obs_YY];
obs_2d = sortrows([obs_2d, theta], 3); obs_2d = obs_2d(:, 1:2);
obs_2d = [obs_2d; obs_2d(1, :)];

% Interpolate to grid
interpolant_Ux = scatteredInterpolant(XX,YY,Ux,'natural','none');
interpolant_Uy = scatteredInterpolant(XX,YY,Uy,'natural','none');

% Grid
[xx,yy] = meshgrid(linspace(-1e-3,1e-3,1001), linspace(-4e-4,4e-4,1001)); 
% The 501th row/column is position ZERO.
profile_x_sim = linspace(-1e-3,1e-3,1001) * 1e6; % unit: um
profile_y_sim = linspace(-4e-4,4e-4,1001) * 1e6;

% Interpolate
Ux_interp = interpolant_Ux(xx,yy);
Uy_interp = interpolant_Uy(xx,yy);

% calculate the U_max,ch
[U_max_xx,U_max_yy] = meshgrid(linspace(-7e-4,-3e-4,100), 0);  
U_max_ch = mean(interpolant_Ux(U_max_xx,U_max_yy));
profile_ux_x_sim = Ux_interp(501, :)/U_max_ch;
profile_ux_y_sim = Ux_interp(:, 501)/U_max_ch;

% % use alpha shape for mask
% shp = alphaShape(XX,YY,2e-3/800);
% figure
% plot(shp)
% Ux_interp(~inShape(shp,xx,yy)) = NaN;

figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);

t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
nexttile

contourf(xx,yy,Ux_interp/U_max_ch,100,'LineStyle','none'); hold on;
shading interp
axis equal; axis off

xlim([-4 4]*1e-4)
ylim([-4 4]*1e-4)

colormap('jet')
caxis([0 1.5])

% Plot streanlines
hold on
startX = -4e-4 * ones(1, 50);
startY = linspace(-4e-4,4e-4,50);
lineobj = streamline(xx,yy,Ux_interp/U_max_ch,Uy_interp/U_max_ch, startX, startY);
for ii = 1:numel(startX)
    lineobj(ii).LineWidth = 0.2;
    lineobj(ii).Color = 'k';
end

hold on
plot([0 0], [-7e-4,7e-4], 'Color', [228,26,28]/255, 'LineWidth', 2); hold on
plot([-7e-4,7e-4], [0 0], 'Color', [228,26,28]/255, 'LineWidth', 2); hold on
text(-4e-4,4.3e-4,'Simulation','FontSize',18)

% load uPIV data (x, y are switched compared to the simulation)
load(['F:\Processing & Results\PIV & PTV\20221031-uPIV\ave_V\Round2_Lens=10_' ...
    'Tracer=1.1um_FlowRate=11nL_Atte=50_Dt=8000us_Z=24.08um_aveV.mat']);
mask = ~logical(ave_field_phy.vy); 
mask = xor(bwareafilt(mask, 3), bwareafilt(mask, 2)); % to select the 3rd big pattern in the mask (the obstacle)
obs_center = round(regionprops(mask).Centroid); 
obs_center(2) = obs_center(2) - 1; % correct the obstacle position

profile_x_exp = ave_field_phy.y - ave_field_phy.y(obs_center(1));
profile_y_exp = ave_field_phy.x - ave_field_phy.x(obs_center(2)); % make the center position be zero
real_channel_halfWidth = 390; expect_channel_halfWidth = 400; % these are estimations
profile_y_exp = profile_y_exp / real_channel_halfWidth  * expect_channel_halfWidth; % calibration

profile_ux_x_exp = ave_field_phy.vy(obs_center(2), :);
profile_ux_y_exp = ave_field_phy.vy(:, obs_center(1));
U_max_ch_exp = max(profile_ux_x_exp);
profile_ux_x_exp = profile_ux_x_exp/U_max_ch_exp;
profile_ux_y_exp = profile_ux_y_exp/U_max_ch_exp;


nexttile

contourf(ave_field_phy.y, ave_field_phy.x, ave_field_phy.vy/U_max_ch_exp, 100,'LineStyle','none'); 
hold on; axis equal;  axis off

xlim([ave_field_phy.y(obs_center(1))-real_channel_halfWidth, ave_field_phy.y(obs_center(1))+real_channel_halfWidth])
ylim([ave_field_phy.x(obs_center(2))-real_channel_halfWidth, ave_field_phy.x(obs_center(2))+real_channel_halfWidth])

colormap('jet')
caxis([0 1.5])

c = colorbar;
c.Label.String = '$u_x/U_{\rm max,ch}$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 18;

% Plot streanlines
hold on
startX = ave_field_phy.y(obs_center(1))-real_channel_halfWidth * ones(1, 50);
startY = linspace(ave_field_phy.x(obs_center(2))-real_channel_halfWidth, ave_field_phy.x(obs_center(2))+real_channel_halfWidth,50);
lineobj = streamline(ave_field_phy.y,ave_field_phy.x,ave_field_phy.vy/U_max_ch_exp, ave_field_phy.vx/U_max_ch_exp, startX, startY);
for ii = 1:numel(startX)
    lineobj(ii).LineWidth = 0.2;
    lineobj(ii).Color = 'k';
end

hold on
line([ave_field_phy.y(obs_center(1)) ave_field_phy.y(obs_center(1))], [0 900], 'Color', [77,175,74]/255, 'LineWidth', 2); hold on
line([0 900], [ave_field_phy.x(obs_center(2)) ave_field_phy.x(obs_center(2))], 'Color', [77,175,74]/255, 'LineWidth', 2); hold on
text(35,945,'Experiment','FontSize',18)

f=gcf;
savefig(f,'D:\Dropbox\Research\My PhD thesis\Figures\2-methods\PIV\uPIV20221031_OpenFOAM_TriObs_flowfield.fig')
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-painters','D:\Dropbox\Research\My PhD thesis\Figures\2-methods\PIV\uPIV20221031_OpenFOAM_TriObs_flowfield.eps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare the velocity profiles
figure('color', 'w');  set_plot(gcf, gca)
plot(profile_x_exp, profile_ux_x_exp,'Color', [77,175,74]/255, 'LineWidth', 2, ...
    'LineStyle','none', 'Marker','o','MarkerSize', 6); hold on
plot(profile_x_sim, profile_ux_x_sim, 'color', [228,26,28]/255, 'LineWidth', 2);
xlim([-400 400]); ylim([0 1.5])
xlabel('$x\ (\mathrm{\mu m})$');
ylabel('$u_x/U_{\rm max, ch}$');
legend('$\rm{Experiment}$','$\rm{Simulation}$','location','northeast')

f=gcf;
savefig(f,'D:\Dropbox\Research\My PhD thesis\Figures\2-methods\PIV\uPIV20221031_OpenFOAM_TriObs_Ux_X.fig')
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-painters','D:\Dropbox\Research\My PhD thesis\Figures\2-methods\PIV\uPIV20221031_OpenFOAM_TriObs_Ux_X.eps')

figure('color', 'w');  set_plot(gcf, gca)
plot(profile_y_exp+expect_channel_halfWidth, profile_ux_y_exp,'Color', [77,175,74]/255, 'LineWidth', 2, ...
    'LineStyle','none', 'Marker','o','MarkerSize', 6); hold on
plot(profile_y_sim+expect_channel_halfWidth, profile_ux_y_sim, 'color', [228,26,28]/255, 'LineWidth', 2);
xlim([0 800]); ylim([0 2])
xlabel('$y\ (\mathrm{\mu m})$');
ylabel('$u_x/U_{\rm max, ch}$');

f=gcf;
savefig(f,'D:\Dropbox\Research\My PhD thesis\Figures\2-methods\PIV\uPIV20221031_OpenFOAM_TriObs_Ux_Y.fig')
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-painters','D:\Dropbox\Research\My PhD thesis\Figures\2-methods\PIV\uPIV20221031_OpenFOAM_TriObs_Ux_Y.eps')