%% Select the 2nd frame in the uPIV (the one contains the fiber)
clear; close all; clc

set(0, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

selpath = uigetdir(['Z:\Experimental Data (RAW)\FSI - Rigid Fiber &  Individual ' ...
    'Obstacle\20220817_Channel2022-07-26_UpTrianglePillar_SU8fibre-and-Tracer\']);
listing = dir(selpath);
pathname = uigetdir(['F:\Experimental Data (EXTRACTED)\FSI - Rigid Fiber &  ' ...
    'Individual Obstacle\20220817_Channel2022-07-26_UpTrianglePillar_SU8fibre-and-Tracer\'], 'Choose a folder to save the results');

Names = {listing.name};  % the file/folder names
cases_ind = contains(Names,'Lateral'); % only select cases.
thefolder_ind = cell2mat({listing.isdir}); % only select the folders.
to_be_processed = and(cases_ind, thefolder_ind);
to_be_processed_ind = find(to_be_processed==1); % the indices of the PAs folders.

caseNo = sum(to_be_processed);
for ii = 1:caseNo

    theONE = Names{to_be_processed_ind(ii)};
    d = dir(fullfile(selpath,theONE,'\CreateTimeSeries\*.im7')); nn = length(d); % number of *.im7 files.
    v = loadvec(fullfile(selpath,theONE,['\CreateTimeSeries\B[2:2:',num2str(nn),'].im7']));  % to calculate the ave_field
    v_subaveBG = subaverf(v);  % subtract the mean field

    newimg = zeros(length(v(1).x), length(v(1).y), nn/2);
    newimg = uint16(newimg);

    if ~exist([pathname,filesep,'AfterAveBGR'],'dir')
        mkdir([pathname,filesep,'AfterAveBGR'])
    end

    for jj = 1:nn/2
        newimg(:, :, jj) = uint16((v_subaveBG(jj).w)*2^4);
    end

    options.append = true;
    saveastiff(newimg, [pathname,filesep,'AfterAveBGR',filesep,theONE,'_AABGR.tif'], options);

end