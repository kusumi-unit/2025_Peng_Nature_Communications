%% Initial setteings
% Calculate the KD value from the PCCF 
%
% G + G <--> GG or A + A <--> AA: kon1, koff1
% R + R <--> RR or B + B <--> BB: kon2, koff2
% G + R <--> GR or A + B <--> AB: kon3, koff3

clear;

% Pre-determined parameters from homo-dimerization experiment
kon1 = 1.143625; 
koff1 = 6.711409;
kon2 = 0.487256;
koff2 = 8.000000; 
% Pre-determined parameters from hetero-dimerization lifetime experiment
koff3 = 4.149377593;

% Results for the Gaussian fitting of PCCF 
% y = y0 + (C/(w*sqrt(pi/2)))*exp(-2*((x-xc)/w)^2)
C = 0.95635; 
W = 85.62388;
C_SE = 0.03106; % SE for parameter
W_SE = 3.21136; % SE for parameter

% Experimetal conditions
Area = 100; % observation area, um^2
NtRealG_M = 71.1725; % # of observed green spots / frame
NtRealR_M = 73.223; % # of observed red spots / frame
EG = 0.663087205; % Labeling efficiency for green spots: 0-1
ER = 0.858544064; % Labeling efficiency for green spots: 0-1
InciTh = 200; % threshold distance for incidental colocalizations (nm)

% Dummy input
Inci = 1;
Homo = 1;
DobsG = 0;
DobsR = 0;


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
    ParaEst2_Na_b_KDcalc(EstNab,NtRealG_M,NtRealR_M,DobsG,DobsR,kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER,Area);


%% Error estimation
% Upper error limit calculation
C_U = C - C_SE;
W_U = W - W_SE;

if Inci == 1
    F3_U = -(C_U*NtReal^2*pi*W_U^2)/(2*A);
end

EstNab_U = (-F2+sqrt(F2^2-4*F1*F3_U))/(2*F1);
[EstNa_U,EstNb_U,EstNaa_U,EstNbb_U,Estkon3_U,EstKd_U]=...
    ParaEst2_Na_b_KDcalc(EstNab_U,NtRealG_M,NtRealR_M,DobsG,DobsR,kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER,Area);

% Lower error limit calculation
C_L = C + C_SE;
W_L = W + W_SE;

if Inci == 1
    F3_L = -(C_L*NtReal^2*pi*W_L^2)/(2*A);
end

EstNab_L = (-F2+sqrt(F2^2-4*F1*F3_L))/(2*F1);
[EstNa_L,EstNb_L,EstNaa_L,EstNbb_L,Estkon3_L,EstKd_L]=...
    ParaEst2_Na_b_KDcalc(EstNab_L,NtRealG_M,NtRealR_M,DobsG,DobsR,kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER,Area);

% Sum error
EstKd_SE = (EstKd_U - EstKd_L)/2;

% Summarize the table
%[KD3, KD3_SE, Kon3, Kon3_SE, Na, Na_SE, Nb, Nb_SE, Naa, Naa_SE Nbb, Nbb_SE, Nab, Nab_SE, TG, TR]
Sum_results(1,1) = EstKd;
Sum_results(1,2) = EstKd_SE;
Sum_results(1,3) = Estkon3;
Sum_results(1,4) = abs(Estkon3_U - Estkon3_L)/2;
Sum_results(1,5) = EstNa;
Sum_results(1,6) = abs(EstNa_U - EstNa_L)/2;
Sum_results(1,7) = EstNb;
Sum_results(1,8) = abs(EstNb_U - EstNb_L)/2;
Sum_results(1,9) = EstNaa;
Sum_results(1,10) = abs(EstNaa_U - EstNaa_L)/2;
Sum_results(1,11) = EstNbb;
Sum_results(1,12) = abs(EstNbb_U - EstNbb_L)/2;
Sum_results(1,13) = EstNab;
Sum_results(1,14) = abs(EstNab_U - EstNab_L)/2;
Sum_results(1,15) = NtRealG_M;
Sum_results(1,16) = NtRealR_M;

header(1,1) = "KD";
header(1,2) = "SE";
header(1,3) = "kon3";
header(1,4) = "SE";
header(1,5) = "# of A monomer / frame";
header(1,6) = "SE";
header(1,7) = "# of B monomer / frame";
header(1,8) = "SE";
header(1,9) = "# of AA dimer / frame";
header(1,10) = "SE";
header(1,11) = "# of BB dimer / frame";
header(1,12) = "SE";
header(1,13) = "# of AB dimer / frame";
header(1,14) = "SE";
header(1,15) = "True # of A protomer / frame";
header(1,16) = "True # of B protomer / frame";

Sum_results = vertcat(header,Sum_results);

writematrix(Sum_results,"Results.xlsx");