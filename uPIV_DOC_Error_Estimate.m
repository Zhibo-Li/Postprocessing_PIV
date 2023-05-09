clear; close all; clc

constant=struct('Kb',1.380649*10^(-23),'e',1.602176*10^(-19),...
    'Na',6.02214076*10^23,'epsilon0',8.8541878128*10^(-12),'c',299792458,...
    'ThermalEnergy',4.11*10^(-21),'visc_water_23K',0.00092353,'Plank',6.62607004*10^(-34));

SI=struct('meter',1,'second',1,'minute',60,'hour',3600,'newton',1,'joule',1,'kg',1,'g',10^-3,'liter',10^-3,'hertz',1,...
    'kelvin',1,'watt',1,'coulomb',1,'volt',1,'farad',1,'ohm',1,'pascal',1,'bar',10^5,...
    'angstrom',10^-10,'eV',1.602176565*10^-19,'zepto',10^-21,'atto',10^-18,'femto',10^-15,...
    'pico',10^-12,'nano',10^-9,'micro',10^-6,'milli',10^-3,'centi',10^-2,'deci',10^-1,'kilo',10^3,'mega',10^6,'giga',10^9,...
'percent',10^-2);


% input 1
% d_p = 1.1e-6; % particle diameter, unit: m
% Temp = 300; % temperature, unit: K
% UU = 600e-6; % typical velocity, unit: m/s
% tt = 10000e-6; % typical time, unit: s
% NN = 500*3; % particle number 
% mu = 340e-3; % viscosity, unit: Pa.s
% 
% epsilon_DOC = 0.01;
% NA = 0.25; % numerical aperture
% MM = 10; % magnitude
% n_0 = 1; % refractive index
% lambda_DOC = 532e-9; % wavelength, unit: m

% input 2
d_p = 0.3e-6; % particle diameter, unit: m
Temp = 300; % temperature, unit: K
UU = 10e-6; % typical velocity, unit: m/s
tt = 10000e-6; % typical time, unit: s
NN = 500*5; % particle number 
mu = 6.1e-3; % viscosity, unit: Pa.s

epsilon_DOC = 0.01;
NA = 1.2; % numerical aperture
MM = 63; % magnitude
n_0 = 1.33; % refractive index
lambda_DOC = 532e-9; % wavelength, unit: m


% calculation
error_B = 1/UU * sqrt(2*constant.Kb*Temp/3/pi/mu/d_p/tt);
error_B_ = error_B/sqrt(NN);

DOC = 2 * sqrt((1-sqrt(epsilon_DOC))/sqrt(epsilon_DOC)*(n_0^2*d_p^2/4/NA^2 + ...
    5.95*(MM+1)^2*lambda_DOC^2*n_0^2/16/MM^2/NA^4)) * 1e6; % unit: m