function [EstNab,EstNa,EstNb,EstKd,EstKd_SE,FitResults] = ParaEst(Area,NtRealG_M,NtRealR_M,EG,ER,InciTh,GaussResults,Inci,Homo)

% for debug
%{
clear all
Area = 100;
NtRealG_M = 70.89923077;%86.283;
NtRealR_M = 69.79923077;%85.2215;
EG = 0.74;
ER = 0.7;
C = 17.80857;
W = 68.04102;
C_SE = 0.62101;
W_SE = 2.08773;
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

% Convert the unit
A = Area*10^6;
% Calibrate incidental fraction
Th = InciTh;
if Inci == 1
    Napp_Gall = (A-A*sqrt(1-2*pi*Th^2*NtRealG_M/A))/(pi*Th^2);
    Napp_Rall = (A-A*sqrt(1-2*pi*Th^2*NtRealR_M/A))/(pi*Th^2);
    Napp_Tall = Napp_Gall + Napp_Rall;
    F1 = A/(2*pi*NtRealG_M*NtRealR_M*W^2);
    F2 = Th^2/(4*W^2*NtRealG_M*NtRealR_M);
    Ninc_TG = pi*Th^2*Napp_Gall^2/(2*A);
    Ninc_TR = pi*Th^2*Napp_Rall^2/(2*A);
else
end

% # of dimer estimation
EstNab = (2*F2*Napp_Tall-F1+sqrt((F1-2*F2*Napp_Tall)^2+12*F2*C))/(6*F2*EG*ER);
% # of monomer extimation
EstNa = (NtRealG_M + (EG^2/2-2*EG)*EstNab + Ninc_TG)/EG;
EstNb = (NtRealR_M + (ER^2/2-2*ER)*EstNab + Ninc_TR)/ER;


if Homo == 1
    % Homo-dimer existing case---------------------------------------------
    % Kd estimation
    EstKd = (EstNa+EstNb)^2/(2*EstNab*A)*10^6; % convert nm to um
    % Error estimation
    W2 = W^2;
    W2_SE = sqrt(2*W*W_SE); 
    CW2_SE = sqrt((W2*C_SE)^2+(C*W2_SE)^2);
    Nab_SE = (NtRealG_M*NtRealR_M*2*pi)/(EG*ER*A)*CW2_SE;
    Na_SE = (EG^2/2-2*EG)/(EG)*Nab_SE;
    Nb_SE = (ER^2/2-2*ER)/(ER)*Nab_SE;
    
    Na_b_SE = sqrt(Na_SE^2 + Nb_SE^2);
    Nume = (EstNa+EstNb)^2*10^6;
    Nume_SE = sqrt(2*(EstNa+EstNb)*Na_b_SE)*10^6;
    Deno = (2*EstNab*A);
    Deno_SE = 2*A*Nab_SE;
    
    EstKd_SE = sqrt((Nume_SE/Deno)^2 + (Nume/Deno^2*Deno_SE)^2);
    
else
    % Currently not available
    %{
    % No homo-dimer existing case------------------------------------------
    EstKd = ((NtApp/2)-EstNab)^2/(EstNab*A)*10^6;
    % Error estimation
    %{
    Upper = ((NtApp/2)-EstNab)^2;
    SE_Upper = 2*Upper*SE_Nab;
    Lower = (EstNab*A);
    SE_Lower = SE_Nab*A;
    Kd_SE = sqrt((1/Lower*SE_Upper)^2+(Upper/Lower^2*SE_Lower)^2)*10^6;
    %}
    %}
end

% Re-shape fitting results
FitResults(1) = C; FitResults(2) = C_SE; 
FitResults(3) = W; FitResults(4) = W_SE;
FitResults(5) = GaussResults(3,1); FitResults(6) = GaussResults(3,2);

end