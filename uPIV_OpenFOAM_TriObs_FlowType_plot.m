%% flowType plot
%% OpenFOAM (202302_triangular_pillar_EXPconfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% read simulation data
data = readmatrix(['F:\Simulation\202302_triangular_pillar_EXPconfig\Data\' ...
    'MidPlane_FlowType.csv']);
flowType = data(1:end, 1);
XX = data(1:end, 2);
YY = data(1:end, 3);

% Interpolate to grid
interpolant_flowType = scatteredInterpolant(XX,YY,flowType,'natural','none');

% Grid
Grid_density = 1001;
Grid_x1 = -5e-4; Grid_x2 = 5e-4; Grid_y1 = -4.5e-4; Grid_y2 = 5.5e-4;
[xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

% Interpolate
flowType_interp = interpolant_flowType(xx,yy);

% To plot
figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
% Plot contour: flowType
contourf(xx,yy,flowType_interp,100,'LineStyle','none');
shading interp
axis equal; axis off

% Plot pillars
hold on
plot([-43.3e-6 0], [-25e-6 50e-6], 'k');
plot([43.3e-6 0], [-25e-6 50e-6], 'k');
plot([-43.3e-6 43.3e-6], [-25e-6 -25e-6], 'k');

xlim([-1 1]*1e-4)
ylim([-1 1]*1e-4)

cmocean('balance');
caxis([-1 1])
% c = colorbar;
% c.Label.String = '$flowType$';
% c.Label.Interpreter = 'LaTeX';
% c.TickLabelInterpreter = 'LaTeX';
% c.FontSize = 18;




%% uPIV (Exp: 20221031)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

[filename, pathname] = uigetfile('F:\Processing & Results\PIV & PTV\20221031-uPIV\ave_V\*.mat');
load(fullfile(pathname,filename))

%%%%% calculate flow type indicator %%%%%
u = ave_field_phy.vx; v = ave_field_phy.vy;
x = ave_field_phy.x; y = ave_field_phy.y;
[Y,X] = meshgrid(y, x);

dx = x(2) - x(1); dy = y(2) - y(1); % here dx == dy

[du_dy,du_dx] = gradient(u); % x, u: longer edge;  y, v: shorter edge.
[dv_dy,dv_dx] = gradient(v);

du_dx = du_dx / dx; dv_dx = dv_dx / dx;
du_dy = du_dy / dy; dv_dy = dv_dy / dy;

VGT_sym = [du_dx(:), 0.5*(du_dy(:)+dv_dx(:)), dv_dy(:), 0.5*(du_dy(:)+dv_dx(:))];
VGT_anti = [zeros(numel(du_dx),1), 0.5*(dv_dx(:)-du_dy(:)), zeros(numel(du_dx),1), 0.5*(du_dy(:)-dv_dx(:))];

% flow type indicator (lambda)
flow_type = nan(numel(du_dx),1);
for jj = 1: numel(du_dx)

    norm_E = norm(reshape(VGT_sym(jj, :), 2, 2),'fro');
    norm_O = norm(reshape(VGT_anti(jj, :), 2, 2),'fro');

    flow_type(jj) = (norm_E-norm_O) / (norm_E+norm_O);

end
flow_type = reshape(flow_type, size(X, 1), size(X, 2));


% To plot
figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
% Plot contour: flowType
contourf(X, Y, flow_type,100,'LineStyle','none');
shading interp
axis equal; axis off

xlim([415 632])
ylim([310 550])

cmocean('balance');
caxis([-1 1])




