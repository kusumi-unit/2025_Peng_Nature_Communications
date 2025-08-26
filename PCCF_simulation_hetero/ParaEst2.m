function [EstNab,EstNa,EstNb,EstNaa,EstNbb,Estkon3,EstKd,EstKd_SE,FitResults]...
    = ParaEst2(Area,NtRealG_M,NtRealR_M,DobsG,DobsR,EG,ER,InciTh,GaussResults,...
    Inci,Homo,kon1,koff1,kon2,koff2,koff3)

% G + G <--> GG or A + A <--> AA: kon1, koff1
% R + R <--> RR or B + B <--> BB: kon2, koff2
% G + R <--> GR or A + B <--> AB: kon3, koff3

% for debug
%{
clear all;
% Kd: 3
% True GR: 11.4616
% True kon3: 2.666
% True kon3: 
% Ntrue G all: 75.30
% Ntrue R all: 75.3
% Ntrue GG: 2.11;
% Ntrue RR: 4.86;

kon1 = 0.47; % homo-DOR before stimlation
koff1 = 8.00;
kon2 = 1.03; % homo-KOR before stimulation
koff2 = 6.71;
koff3 = 8; % hetero-dimer

Area = 100;
NtRealG_M = 45;
NtRealR_M = 45;
EG = 0.74;
ER = 0.87;
C = 2.6153;
W = 133.4508;
C_SE = 0.024216;
W_SE = 2.6124;
InciTh = 200;
Inci = 1;
Homo = 1;
%}

% pickup fitting data
%
C = GaussResults(1,1);
W = GaussResults(2,1);
C_SE = GaussResults(1,2);
W_SE = GaussResults(2,2);
%}


%% Calculation

% Convert the unit
A = Area*10^6;
% Calibrate incidental fraction
Th = InciTh;
if Inci == 1
    Napp_Gall = (A-A*sqrt(1-2*pi*Th^2*NtRealG_M/A))/(pi*Th^2);
    Napp_Rall = (A-A*sqrt(1-2*pi*Th^2*NtRealR_M/A))/(pi*Th^2);
    NtReal = NtRealG_M + NtRealG_M;
    F1 = (3*pi*Th^2*EG^2*ER^2)/(2*A);
    F2 = EG*ER*(1-(pi*Th^2)/A*(Napp_Gall+Napp_Rall));
    F3 = -(C*NtReal^2*pi*W^2)/(2*A);    
    Ninc_TG = pi*Th^2*Napp_Gall^2/(2*A);
    Ninc_TR = pi*Th^2*Napp_Rall^2/(2*A);
else
end

% # of hetero-dimer estimation
EstNab = (-F2+sqrt(F2^2-4*F1*F3))/(2*F1);
% % of monomer, homo-dimer, kon3 estimation; should be done in um unit
[EstNa,EstNb,EstNaa,EstNbb,Estkon3,EstKd]=...
    ParaEst2_Na_b(EstNab,NtRealG_M,NtRealR_M,DobsG,DobsR,kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER,Area);


%% Error estimation
% Upper error limit calculation
C_U = C - C_SE;
W_U = W - W_SE;

if Inci == 1
    F3_U = -(C_U*NtReal^2*pi*W_U^2)/(2*A);
end

EstNab_U = (-F2+sqrt(F2^2-4*F1*F3_U))/(2*F1);
[EstNa_U,EstNb_U,EstNaa_U,EstNbb_U,Estkon3_U,EstKd_U]=...
    ParaEst2_Na_b(EstNab_U,NtRealG_M,NtRealR_M,DobsG,DobsR,kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER,Area);

% Lower error limit calculation
C_L = C + C_SE;
W_L = W + W_SE;

if Inci == 1
    F3_L = -(C_L*NtReal^2*pi*W_L^2)/(2*A);
end

EstNab_L = (-F2+sqrt(F2^2-4*F1*F3_L))/(2*F1);
[EstNa_L,EstNb_L,EstNaa_L,EstNbb_L,Estkon3_L,EstKd_L]=...
    ParaEst2_Na_b(EstNab_L,NtRealG_M,NtRealR_M,DobsG,DobsR,kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER,Area);

% Sum error
EstKd_SE = (EstKd_U - EstKd_L)/2;

% If the error can't estimated, put nan
if isnan(EstKd_SE)
    EstNa=nan;EstNb=nan;EstNaa=nan;EstNbb=nan;Estkon3=nan;EstKd=nan;
end

% Re-shape fitting results
FitResults(1) = C; FitResults(2) = C_SE; 
FitResults(3) = W; FitResults(4) = W_SE;
FitResults(5) = GaussResults(3,1); FitResults(6) = GaussResults(3,2);

end