%% Pillar array area
clear; close all; clc

set(0, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

selpath = uigetdir('G:\PhD, PMMH, ESPCI\Experimental Data (EXTRACTED)\');
listing = dir(selpath);

Names = {listing.name};  % the file/folder names
XblankArea = xor(contains(Names,'Lens='), contains(Names,'BlankArea_')); % only select the PAs related files/folders.
thefolder_ind = cell2mat({listing.isdir}); % only select the folders.
PAsArea = and(XblankArea, thefolder_ind); % the PAs folders.
PAsArea_ind = find(PAsArea==1); % the indices of the PAs folders.

magni = 0.067; % the magnification of the objective (um/pixel)

% get the cross-points among four adjacent pillars (in pixel). [Manually input]
tiltshift = 2; thegap = 632;
x1 = 491; y1 = 318;    x2 = 810; y2 = 630;
theX1 = zeros(3,3); theY1 = zeros(3,3);   theX2 = zeros(3,3); theY2 = zeros(3,3);

theX1(:, 1) = transpose(x1-tiltshift*(0:2)); for jj = 2:3; theX1(:, jj) = theX1(:, jj-1) + thegap; end 
theX2(:, 1) = transpose(x2-tiltshift*(0:2)); for jj = 2:3; theX2(:, jj) = theX2(:, jj-1) + thegap; end 

theY1(1, :) = y1+tiltshift*(0:2); for jj = 2:3; theY1(jj, :) = theY1(jj-1, :) + thegap; end 
theY2(1, :) = y2+tiltshift*(0:2); for jj = 2:3; theY2(jj, :) = theY2(jj-1, :) + thegap; end 

theCrossPoints = [reshape(complex(theX1, theY1), 1, []), reshape(complex(theX2, theY2), 1, [])];

% normalize the crosspoints by the calculation gap length (for V-profile along Z).
LocSelectX = round(real(theCrossPoints)/12); % 12 is the calculation gap length.
LocSelectY = round(imag(theCrossPoints)/12);
crossnum = length(theCrossPoints);

% normalize the X by the calculation gap length (for V-profile along flow direction).
theNormX = round([theX1(2, :), theX2(2, :)]/12);

layersNUM = sum(PAsArea);
Z = zeros(1, layersNUM);  UUU_alongZ = zeros(1, layersNUM); VVV_alongZ = zeros(1, layersNUM);
UUU_alongZ_err = zeros(1, layersNUM); VVV_alongZ_err = zeros(1, layersNUM); 
for ii = 1:layersNUM

    theONE = Names{PAsArea_ind(ii)};
    Z(ii) = str2double(extractBetween(theONE,'Z=','um'));
    deltaT = str2double(extractBetween(theONE,'Dt=','us'));
%     v = loadvec(fullfile(selpath,theONE,'\AddGeometricMask_01\PIV_MPd(2x48x48_75%ov)\*.vc7'));  % to calculate the ave_field
%     ave_field = averf(v);
%     [ave_field_AF, ave_field_STD, ave_field_RMS] = averf(v);
    load(fullfile(selpath,theONE,'\AddGeometricMask_01\PIV_MPd(2x48x48_75%ov)\ave_field_STD_RMS.mat')); % load the ave_field
%     save(fullfile(selpath,theONE,'\AddGeometricMask_01\PIV_MPd(2x48x48_75%ov)\ave_field_STD_RMS.mat'), "ave_field_AF", "ave_field_RMS", "ave_field_STD");
    ave_field_phy = ave_field_AF; % in physical units.
    ave_field_phy.vx = ave_field_AF.vx * magni / deltaT; % velocity in m/s.
    ave_field_phy.vy = ave_field_AF.vy * magni / deltaT; % velocity in m/s.
    ave_field_phy.x = ave_field_AF.x * magni; % velocity in um.
    ave_field_phy.y = ave_field_AF.y * magni; % velocity in um.

    ave_field_phy.unitvx = 'm/s'; ave_field_phy.unitvy = 'm/s';
    ave_field_phy.unitx = '$\mu{m}$'; ave_field_phy.unity = '$\mu{m}$';

    % the U and V along the z-direction.
    UU = []; VV = []; % calculate the speed near the cross-points
    for pp = 1:crossnum
        UU = [UU, mean(mean(ave_field_phy.vx(LocSelectX(pp)-1:LocSelectX(pp)+1,...
            LocSelectY(pp)-1:LocSelectY(pp)+1),'omitnan'))];
        VV = [VV, mean(mean(ave_field_phy.vy(LocSelectX(pp)-1:LocSelectX(pp)+1,...
            LocSelectY(pp)-1:LocSelectY(pp)+1),'omitnan'))]; 
    end
    UUU_alongZ(ii) = mean(UU);  VVV_alongZ(ii) = mean(VV); % the averaged the U, V over 18 cross-points.
    UUU_alongZ_err(ii) = std(UU);  VVV_alongZ_err(ii) = std(VV); % the averaged the U, V over 18 cross-points.

    % the U and V along the flow-direction.
    UU = []; VV = [];
    if 25<Z(ii) && Z(ii)<26   % Z = 25.07: middle plane.
        for pp = 1:length(theNormX)
            tmpU = ave_field_phy.vx(theNormX(pp), :);
            tmpV = ave_field_phy.vy(theNormX(pp), :);
%             foo = find(tmpV(2:end-1) ~= 0); ind(pp) = foo(1);  % try to find the shifting length.
            if pp > 3
                shiftlength = -27; % because it's staggered, should shift the array before the average.
                tmpU = circshift(tmpU,shiftlength);
                tmpV = circshift(tmpV,shiftlength);
            end
            UU = [UU; tmpU];
            VV = [VV; tmpV];
        end
    UUU_alongX = mean(UU, 'omitnan');  VVV_alongX = mean(VV, 'omitnan'); % the averaged the U, V over 18 cross-points.
    UUU_alongX_err = std(UU, 'omitnan');  VVV_alongX_err = std(VV, 'omitnan'); % the averaged the U, V over 18 cross-points.
    end

%     figure; showf(ave_field_phy, 'title', ['Z = ',num2str(Z(ii)), ' $\mu{m}$']);
%     ax = gca; ax.Title.Interpreter = 'latex';

end

%%% Velocity profile along the z-direction (absolute value).
% figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
% errorbar(Z, VVV_alongZ*1e6, VVV_alongZ_err*1e6,'LineStyle','none', 'Marker','*','MarkerSize', 5); hold on
% % % the simulation results: 
% data = readmatrix('D:\Dropbox\PROCESS remotely\Simulation\202208\uPIVexp_H52\Ux_Z-direction.csv');
% Speed_Simu = data(1:10:end, 1); ZZ = data(1:10:end, 9);
% plot(ZZ*1e5, Speed_Simu*1e7, 'r--.', 'LineWidth', 2); hold on
% xlabel('$Z\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
% ylabel('$U_x\ ({\mu}m/s)$','FontSize',22,'Interpreter', 'latex');
% % % title('Speed in the Middle Plane','FontSize',48);
% set(gca,'FontSize',20)
% legend('PIV','Simulation','location','best','FontSize',20)
% xlim([-0.5 52.5])
% % f=gcf;
% % exportgraphics(f,'D:\Dropbox\Research\All Plottings\General plots\20220811_PIV-Simulation_PAsArea_H52_Ux-Z.png','Resolution',100)
% 
% 
% %%% Velocity profile along the z-direction (normalized).
% figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
% errorbar(Z/52-0.02, VVV_alongZ*1e6/147, VVV_alongZ_err*1e6/147,'LineStyle','none', 'Marker','*','MarkerSize', 5); hold on
% % % the simulation results: 
% data = readmatrix('D:\Dropbox\PROCESS remotely\Simulation\202208\uPIVexp_H52\Ux_Z-direction.csv');
% Speed_Simu = data(1:10:end, 1); ZZ = data(1:10:end, 9);
% plot(ZZ*1e5/52, Speed_Simu*1e7/157, 'r--.', 'LineWidth', 2); hold on
% xlabel('$Z/H$','FontSize',22,'Interpreter', 'latex');
% ylabel('$U_x/U_{peak}\ in\ blank\ area$','FontSize',22,'Interpreter', 'latex');
% % % title('Speed in the Middle Plane','FontSize',48);
% set(gca,'FontSize',20)
% legend('PIV','Simulation','location','best','FontSize',20)
% xlim([0 1])
% % f=gcf;
% % exportgraphics(f,'D:\Dropbox\Research\All Plottings\General plots\20220811_PIV-Simulation_PAsArea_H52_Ux-Z_normalized.png','Resolution',100)    


%%% Velocity profile along the flow-direction (absolute value).
% XXX_alongX = ave_field_phy.y;
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % the simulation results: 
% data = readmatrix('D:\Dropbox\PROCESS remotely\Simulation\202208\uPIVexp_H52\Ux_X-direction.csv');
% Speed_Simu = data(1:2:end, 1); XX = data(1:2:end, 7);
% errorbar(XXX_alongX(1:140), VVV_alongX(1:140)*1e6, VVV_alongX_err(1:140)*1e6,'LineStyle','none', 'Marker','*','MarkerSize', 5); hold on
% plot(XX*1e5, Speed_Simu*1e7, 'r--.', 'LineWidth', 2); hold on
% xlabel('$X\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
% ylabel('$U_x\ ({\mu}m/s)$','FontSize',22,'Interpreter', 'latex');
% % % title('Speed in the Middle Plane','FontSize',48);
% set(gca,'FontSize',20)
% legend('PIV','Simulation','location','best','FontSize',20)
% xlim([-0.5 115])
% % f=gcf;
% % exportgraphics(f,'D:\Dropbox\Research\All Plottings\General plots\20220811_PIV-Simulation_PAsArea_H52_Ux-X.png','Resolution',100)   


%%% Velocity profile along the flow-direction (normalized).
% XXX_alongX = ave_field_phy.y;
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % the simulation results: 
% data = readmatrix('D:\Dropbox\PROCESS remotely\Simulation\202208\uPIVexp_H52\Ux_X-direction.csv');
% Speed_Simu = data(1:2:end, 1); XX = data(1:2:end, 7);
% errorbar(XXX_alongX(1:140), VVV_alongX(1:140)*1e6/147, VVV_alongX_err(1:140)*1e6/147,'LineStyle','none', 'Marker','*','MarkerSize', 5); hold on
% plot(XX*1e5, Speed_Simu*1e7/157, 'r--.', 'LineWidth', 2); hold on
% xlabel('$X\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
% ylabel('$U_x/U_{peak}\ in\ blank\ area$','FontSize',22,'Interpreter', 'latex');
% % % title('Speed in the Middle Plane','FontSize',48);
% set(gca,'FontSize',20)
% legend('PIV','Simulation','location','best','FontSize',20)
% xlim([-0.5 115])
% % f=gcf;
% % exportgraphics(f,'D:\Dropbox\Research\All Plottings\General plots\20220811_PIV-Simulation_PAsArea_H52_Ux-X_normalized.png','Resolution',100)  

%% Blank area
clear; close all; clc

selpath = uigetdir('G:\PhD, PMMH, ESPCI\Experimental Data (EXTRACTED)\');
listing = dir(selpath);

Names = {listing.name};  % the file/folder names
blankArea = and(contains(Names,'Lens='), contains(Names,'BlankArea_')); % only select the Blank-area related files/folders.
thefolder_ind = cell2mat({listing.isdir}); % only select the folders.
PAsArea = and(blankArea, thefolder_ind); % the Blank-area folders.
PAsArea_ind = find(PAsArea==1); % the indices of the Blank-area folders.

magni = 0.067; % the magnification of the objective (um/pixel)

layersNUM = sum(PAsArea);
Z = zeros(1, layersNUM);  UUU_alongZ = zeros(1, layersNUM); VVV_alongZ = zeros(1, layersNUM);
UUU_alongZ_err = zeros(1, layersNUM); VVV_alongZ_err = zeros(1, layersNUM); 
for ii = 1:layersNUM

    theONE = Names{PAsArea_ind(ii)};
    Z(ii) = str2double(extractBetween(theONE,'Z=','um'));
    deltaT = str2double(extractBetween(theONE,'Dt=','us'));
%     v = loadvec(fullfile(selpath,theONE,'\AddGeometricMask\PIV_MPd(2x48x48_75%ov)\*.vc7'));  % to calculate the ave_field
%     ave_field = averf(v);
%     [ave_field_AF, ave_field_STD, ave_field_RMS] = averf(v);
    load(fullfile(selpath,theONE,'\AddGeometricMask\PIV_MPd(2x48x48_75%ov)\ave_field_STD_RMS.mat')); % load the ave_field
%     save(fullfile(selpath,theONE,'\AddGeometricMask\PIV_MPd(2x48x48_75%ov)\ave_field_STD_RMS.mat'), "ave_field_AF", "ave_field_RMS", "ave_field_STD");
    ave_field_phy = ave_field_AF; % in physical units.
    ave_field_phy.vx = ave_field_AF.vx * magni / deltaT; % velocity in m/s.
    ave_field_phy.vy = ave_field_AF.vy * magni / deltaT; % velocity in m/s.
    ave_field_phy.x = ave_field_AF.x * magni; % velocity in um.
    ave_field_phy.y = ave_field_AF.y * magni; % velocity in um.

    ave_field_phy.unitvx = 'm/s'; ave_field_phy.unitvy = 'm/s';
    ave_field_phy.unitx = '$\mu{m}$'; ave_field_phy.unity = '$\mu{m}$';

    % the U and V along the y-direction in the middle plane.
    if ii == 1
        counting = 0; 
        UUU_alongY = zeros(1, size(ave_field_phy.vx, 1));
        VVV_alongY = zeros(1, size(ave_field_phy.vx, 1));
        UUU_alongY_err = zeros(1, size(ave_field_phy.vx, 1));
        VVV_alongY_err = zeros(1, size(ave_field_phy.vx, 1));
    end
    if 20<Z(ii) && Z(ii)<30
        UUU_alongY = UUU_alongY + mean(ave_field_phy.vx(:, 5: end-5), 2)'; 
        VVV_alongY = VVV_alongY + mean(ave_field_phy.vy(:, 5: end-5), 2)';
        UUU_alongY_err = UUU_alongY + std(ave_field_phy.vx(:, 5: end-5), 0, 2)'; 
        VVV_alongY_err = VVV_alongY + std(ave_field_phy.vy(:, 5: end-5), 0, 2)';
        counting = counting + 1;
    end
        

%     % the U and V along the z-direction.
    select_U = ave_field_phy.vx(192:201, :);
    select_V = ave_field_phy.vy(192:201, :);
    UUU_alongZ(ii) = mean(mean(select_U));
    VVV_alongZ(ii) = mean(mean(select_V)); 
    UUU_alongZ_err(ii) = std(select_U(:));
    VVV_alongZ_err(ii) = std(select_V(:));

%     figure; showf(ave_field_phy, 'title', ['Z = ',num2str(Z(ii)), ' $\mu{m}$']);
%     ax = gca; ax.Title.Interpreter = 'latex';
end

%%% Velocity profile along the z-direction (absolute value).
% figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
% errorbar(Z, VVV_alongZ*1e6,VVV_alongZ_err*1e6, 'LineStyle','none', 'Marker','*','MarkerSize', 5); hold on
% % % % the simulation results: 
% data = readmatrix('D:\Dropbox\PROCESS remotely\Simulation\202208\uPIVexp_H52\Ux_Z-direction_BlankArea.csv');
% Speed_Simu = data(1:10:end, 1); ZZ = data(1:10:end, 9);
% plot(ZZ*1e5, Speed_Simu*1e7, 'r--.', 'LineWidth', 2); hold on
% % % % theory
% a = 26;   % height/2 (um)
% b = 200; % width/2 (um)
% Q = 2 * 10^(6); %flow rate (nL! the number before *)
% y = 0;
% z = -a:(2*a/100):a;
% Sigma1 = 0; Sigma2 = 0;
% for i = 1:2:10^2
%     Sigma1 = Sigma1 + (-1)^((i-1)/2)*(1-(cosh(i*pi*y/(2*a)))/(cosh(i*pi*b/(2*a))))*(cos(i*pi*z/(2*a)))/(i^3);
%     Sigma2 = Sigma2 + (tanh(i*pi*b/(2*a)))/(i^5);
% end
% Foo = 1 / (1 - 192*a*Sigma2/(pi^5*b));
% f = (12*Q)/((pi)^3*a*b)*Sigma1*Foo;
% plot(z+26, f, 'LineWidth', 2, 'Color','y'); hold on
% xlabel('$Z\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
% ylabel('$U_x\ ({\mu}m/s)$','FontSize',22,'Interpreter', 'latex');
% set(gca,'FontSize',20)
% xlim([-0.5 52.5])
% legend('PIV','Simulation','Theory','location','best','FontSize',20)
% % f=gcf;
% % exportgraphics(f,'D:\Dropbox\Research\All Plottings\General plots\20220811_PIV-Simulation-Theory_BlankArea_H52_Ux-Z.png','Resolution',100)


%%% Velocity profile along the spanwise-direction (absolute value).
% figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
% UUU_alongY = UUU_alongY/counting; VVV_alongY = VVV_alongY/counting;
% UUU_alongY_err = UUU_alongY_err/counting; VVV_alongY_err = VVV_alongY_err/counting;
% XXX_alongY = ave_field_phy.x;
% % errorbar(XXX_alongY, VVV_alongY*1e6, VVV_alongY_err*1e6,'LineStyle','none', 'Marker','*','MarkerSize', 5)
% plot(XXX_alongY(1:end-2), VVV_alongY(1:end-2)*1e6,'LineStyle','none', 'Marker','*','MarkerSize', 5); hold on % (end - 2): remove the 0 value.
% % % % theory
% a = 200;   % height/2 (um)
% b = 26; % width/2 (um)
% Q = 2 * 10^(6); %flow rate (nL! the number before *)
% y = 0;
% z = -a:(2*a/100):a;
% Sigma1 = 0; Sigma2 = 0;
% for i = 1:2:10^2
%     Sigma1 = Sigma1 + (-1)^((i-1)/2)*(1-(cosh(i*pi*y/(2*a)))/(cosh(i*pi*b/(2*a))))*(cos(i*pi*z/(2*a)))/(i^3);
%     Sigma2 = Sigma2 + (tanh(i*pi*b/(2*a)))/(i^5);
% end
% Foo = 1 / (1 - 192*a*Sigma2/(pi^5*b));
% f = (12*Q)/((pi)^3*a*b)*Sigma1*Foo;
% plot(z+212.5, f, 'LineWidth', 2, 'Color','y'); hold on
% xlabel('$Y\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
% ylabel('$U_x\ ({\mu}m/s)$','FontSize',22,'Interpreter', 'latex');
% set(gca,'FontSize',20)
% xlim([-0.5 420.5])
% legend('PIV','Theory','location','best','FontSize',20)
% % f=gcf;
% % exportgraphics(f,'D:\Dropbox\Research\All Plottings\General plots\20220811_PIV-Theory_BlankArea_H52_Ux-Y.png','Resolution',100)



