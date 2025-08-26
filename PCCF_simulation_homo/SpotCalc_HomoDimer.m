function [Dtrue_ALL,Dtrue_GR,Dapp_GR,Dtrue_G,Dtrue_R,Dapp_G,Dapp_R]=SpotCalc_HomoDimer(DtApp,Kd,InciTh,EG,ER,Area)

%{
% for debug
clear all;
DtApp = 2; % Apparent total (Red + Green) spots density (copies / um^2)
Kd = 10; % copies / um^2, 0 for stable dimer, inf for random control
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

% Reverse calculation of the labeling efficiency
E = (EG+ER)/2;
F1 = A*E^2/4 - A*E;
F2 = Napp_ALL + A*E^2*Kd/16;
F3 = 2*F1*F2 - (8*Kd*A^2*E^4)/(16^2);
F4 = F2^2 - (A^2*E^4*Kd^2)/(16^2);

Dtrue_ALL = (-F3 - sqrt(F3^2-4*F1^2*F4))/(2*F1^2);


%% Calculate apprent spot# including labeling efficiency
Dtrue_GR = (Kd + 4*Dtrue_ALL - sqrt(Kd^2+8*Dtrue_ALL*Kd))/16;
Dtrue_GG = Dtrue_GR/2;
Dtrue_RR = Dtrue_GG;

Dtrue_G = (-Kd + sqrt(Kd^2+8*Dtrue_ALL*Kd))/8;
Dtrue_R = Dtrue_G;

Dapp_GR = EG*ER*Dtrue_GR;

Dapp_G = EG*Dtrue_G ...
        + EG^2*Dtrue_GG ...
        + 2*EG*(1-EG)*Dtrue_GG ...
        + EG*(1-ER)*Dtrue_GR;

Dapp_R = ER*Dtrue_R ...
        + ER^2*Dtrue_RR ...
        + 2*ER*(1-ER)*Dtrue_RR ...
        + ER*(1-EG)*Dtrue_GR;

    
end