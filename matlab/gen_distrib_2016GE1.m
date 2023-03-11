close all
clear all

% Add path to functions and load the Granvik model
addpath(genpath('functions'));
%addpath(genpath('period_distribution'));
granvikModel = load("granvikModel/gmb_model.dat");

% Name of the asteroid
name = '2016GE1';

% Rotation period
P        = 0.009438;
P_sigma  = 0.009438*30/100;
% Measurement of the A2 parameter from JPL, units [AU/d^2]
A2       = -1.438865589025615E-12;
sigma_A2 = 4.378E-13;
% Orbital elements
a        = 2.062846039108595;
e        = .5204690304323265;
inc      = 10.72885860239125;
% Convert A2 to da/dt
dadt       = a22my(A2, a, e)
dadt_sigma = a22my(sigma_A2, a, e)
% Absolute magnitude
H       = 26.7;
H_sigma = 0.3;
% Data about diameter
D = nan;
D_sigma = nan;

% Create the folder in the input folder
folder = name; 
cmd = strcat('mkdir ../input/', folder);
eval(cmd);

% Source region probability from Granvik et. al. 2018
[p, s] = gmb_search(a, e, inc, H, granvikModel); 
folder = strcat('input/', folder);
%gen_distrib(D, D_sigma, H, H_sigma, dadt, dadt_sigma, P, P_sigma, p, folder)
gen_distribV2(D, D_sigma, H, H_sigma, dadt, dadt_sigma, P, P_sigma, p, folder)
writeFileAst(name,a,e)
