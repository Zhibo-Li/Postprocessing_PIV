%%% convert *.vc7 to *.mat and convert to physical unit (um/s) and plotting
% direction: x -- spanwise
%            y -- streamwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

magni = 0.42; % the magnification of the objective (um/pixel)

pathname = uigetdir('F:\Processing & Results\PIV & PTV\20221031-uPIV\ave_V\');
% listing = dir(fullfile(pathname, '*.vc7')); suffiex = 'vc7'; % Choose this if you need to convert.
listing = dir(fullfile(pathname, '*.mat')); suffiex = 'mat'; % Choose this if there are *.mat files.

for ii = 1:length(listing)

    theONE = listing(ii).name;
    Z(ii) = str2double(extractBetween(theONE,'Z=','um'));
    deltaT = str2double(extractBetween(theONE,'Dt=','us'));
    if strcmp(suffiex,'vc7')

        v = loadvec(fullfile(pathname,theONE));  % to calculate the physical field

        ave_field_phy = v; % in physical units.
        ave_field_phy.vx = v.vx * magni / deltaT * 1e6; % velocity in um/s.
        ave_field_phy.vy = v.vy * magni / deltaT * 1e6; % velocity in um/s.
        ave_field_phy.x = v.x * magni; % velocity in um.
        ave_field_phy.y = v.y * magni; % velocity in um.
        ave_field_phy.unitvx = '\mu{m}/s'; ave_field_phy.unitvy = '\mu{m}/s';
        ave_field_phy.unitx = '\mu{m}'; ave_field_phy.unity = '\mu{m}';

        save([pathname,filesep,theONE(1:end-4),'.mat'], "ave_field_phy");

    else

        load(fullfile(pathname,theONE))

    end

    %%%%% the U and V along the z-direction (channel height) %%%%%
    if ii == 1
        imshow(imbinarize(ave_field_phy.vy));   
        selrect = getrect; close; % select the area to average.
    end
    select_U = ave_field_phy.vx(selrect(2):selrect(2)+selrect(4), selrect(1):selrect(1)+selrect(3));
    select_V = ave_field_phy.vy(selrect(2):selrect(2)+selrect(4), selrect(1):selrect(1)+selrect(3));
    UUU_alongZ(ii) = mean(mean(select_U));
    VVV_alongZ(ii) = mean(mean(select_V)); 

end

%%% plot V along the z-direction
figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
plot(Z, VVV_alongZ, 'LineStyle','none', 'Marker','*','MarkerSize', 5); hold on
xlabel('$Z\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
ylabel('$U_y\ ({\mu}m/s)$','FontSize',22,'Interpreter', 'latex');
set(gca,'FontSize',20)


%%%%% show U and V in the middle plane %%%%%
[filename, pathname]=uigetfile(fullfile(pathname, '*.mat'), ...
    'Choose the *.mat file of the middle plane');
load(fullfile(pathname, filename));
Z_name = str2double(extractBetween(filename,'Z=','um'));
figure; showf(ave_field_phy, 'title', ['Z = ',num2str(Z_name), ' $\mu{m}$']); % showf function needs PIVMat toolbox.
ax = gca; ax.Title.Interpreter = 'latex';

%%%%% plot V along the y-direction (streamwise) in the middle plane across the obstacle center %%%%%
figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
VV_y = ave_field_phy.vy(round(size(ave_field_phy.vy,1)/2), :);
plot(ave_field_phy.y, VV_y, 'LineStyle','none', 'Marker','*','MarkerSize', 5); 
xlabel('$Y\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
ylabel('$U_y\ ({\mu}m/s)$','FontSize',22,'Interpreter', 'latex');
set(gca,'FontSize',20)

%%%%% plot V along the x-direction (spanwise) in the middle plane across the obstacle center %%%%%
figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
VV_x = ave_field_phy.vy(:, round(size(ave_field_phy.vy,2)/2));
plot(ave_field_phy.x, VV_x, 'LineStyle','none', 'Marker','*','MarkerSize', 5);
xlabel('$X\ ({\mu}m)$','FontSize',22,'Interpreter', 'latex');
ylabel('$U_y\ ({\mu}m/s)$','FontSize',22,'Interpreter', 'latex');
set(gca,'FontSize',20)
