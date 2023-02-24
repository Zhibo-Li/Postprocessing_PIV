clear; close all; clc

pathname = uigetdir('F:\Processing & Results\PIV & PTV\20221031-uPIV\ave_V\');
listing = dir(fullfile(pathname, '*.mat')); suffiex = 'mat'; % Choose this if there are *.mat files.

for ii = 1:length(listing)

    theONE = listing(ii).name;
    load(fullfile(pathname,theONE))

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

    % plot
    pcolor(X, Y, flow_type);shading interp;axis equal;

    % save
    save(fullfile(pathname,theONE), 'ave_field_phy', 'flow_type', 'X', 'Y');

end





