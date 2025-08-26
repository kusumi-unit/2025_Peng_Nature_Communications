function [Dtrue_ALL,Dtrue_GR,Dapp_GR,Dtrue_G,Dtrue_R,Dapp_G,Dapp_R]=SpotCalc_NoDimer(DtApp,InciTh,Area)

%{
% for debug
clear all;
DtApp = 2.5; % Apparent total (Red + Green) spots density (copies / um^2)
Kd = 1; % copies / um^2, 0 for stable dimer, inf for random control
InciTh = 200; % Threshold distance pf incidental homo dimerization (nm)
EG = 0.7; % Labelling efficiancy for Green: 0~1
ER = 0.7;
Area = 100; % um^2
%}


A = Area;
sNT = DtApp * A;
R = InciTh/1000;


%% Calculation for the true total number of spots
% Reverse calculation of the incidental colocalization efect
    % sNT: observed total spots, 
    % corrected by incidental colocalization and labeling efficiency
    % NT_app: apparent spots,
    % corrected by labeling efficiency 
Napp_ALL = 2*(A-A*sqrt(1-pi*R^2*sNT/A))/(pi*R^2);

% No reverse calculation of the labeling efficiency
Dtrue_ALL = Napp_ALL/A;



%% Calculate apprent spot# including labeling efficiency
Dtrue_GR = 0;
Dtrue_GG = 0;
Dtrue_RR = 0;

Dtrue_G = Dtrue_ALL/2;
Dtrue_R = Dtrue_G;

Dapp_GR = 0;

Dapp_G = Dtrue_G;
Dapp_R = Dtrue_R;

    
end
