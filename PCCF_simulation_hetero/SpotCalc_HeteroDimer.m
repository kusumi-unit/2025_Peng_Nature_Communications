function [Dtrue_GALL,Dtrue_RALL,Dtrue_GR,Dtrue_GG,Dtrue_RR,Dapp_GR,Dtrue_G,Dtrue_R,Dapp_G,Dapp_R]=...
    SpotCalc_HeteroDimer(DobsG,DobsR,Kd3,kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER)

% G + G <--> GG or A + A <--> AA: kon1, koff1
% R + R <--> RR or B + B <--> BB: kon2, koff2
% G + R <--> GR or A + B <--> AB: kon3, koff3

%{
% for debug
clear all;
DobsG = 0.711725; % observed total spots density (copies / um^2)
DobsR = 0.73223;
Kd3 = 23.25787522; % hetero-dimer constant (copies / um^2), 0 for stable dimer, inf for random control

%{
kon1 = 0.47; % homo-DOR before stimlation
koff1 = 8.00;
kon2 = 1.03; % homo-KOR before stimulation
koff2 = 6.71;
koff3 = 8; % hetero-dimer
%}
Comb = 'KD'; % MD or KD
if strcmp(Comb,'MD')
    kon1 = 0.53;
    koff1 = 8.47;
    kon2 = 0.47;
    koff2 = 8.00; 
    koff3 = 3.86; 
    FilenameRoot = 'MD_';
elseif strcmp(Comb,'KD')
    kon1 = 1.03; 
    koff1 = 6.71;
    kon2 = 0.47; 
    koff2 = 8.00;
    koff3 = 4.17;
    FilenameRoot = 'KD_';
end
InciTh = 200; % Threshold distance pf incidental homo dimerization (nm)
EG = 0.74; % Labelling efficiancy for Green: 0~1
ER = 0.87;
Area = 100; % um^2
%}


%% Calculation for the true total number of spots
R = InciTh/1000; % convert unit
kon3 = koff3/Kd3;
EA = EG; EB = ER;
% Reverse calculation of the incidental colocalization efect
    % sNT: observed total spots, 
    % corrected by incidental colocalization and labeling efficiency
    % NT_app: apparent spots,
    % corrected by labeling efficiency 
Dapp_A_ALL = (1-sqrt(1-2*pi*R^2*DobsG))/(pi*R^2);
Dapp_B_ALL = (1-sqrt(1-2*pi*R^2*DobsR))/(pi*R^2);


%% Solve simultaneous equations
% Change the name of variables
Da_TA = Dapp_A_ALL;
Da_TB = Dapp_B_ALL;
% Define functions
% TA: true density of toal A 
syms TA TB AA BB AB
eqn1 = kon1*(TA-2*AA-AB)^2-koff1*AA == 0;
eqn2 = kon2*(TB-2*BB-AB)^2-koff2*BB == 0;
eqn3 = kon3*(TA-2*AA-AB)*(TB-2*BB-AB)-koff3*AB == 0;
eqn4 = Da_TA == EA*TA-EA^2*AA;
eqn5 = Da_TB == EB*TB-EB^2*BB;
% Solve functions
S = solve([eqn1,eqn2,eqn3,eqn4,eqn5], [TA,TB,AA,BB,AB]);

% Summarize the result
% V = [TA, TB, AA, BB, AB]
for k=1:length(S.TA)
    V(k,1) = double(vpa(S.TA(k,1)));
end
for k=1:length(S.TB)
    V(k,2) = double(vpa(S.TB(k,1)));
end
for k=1:length(S.AA)
    V(k,3) = double(vpa(S.AA(k,1)));
end
for k=1:length(S.BB)
    V(k,4) = double(vpa(S.BB(k,1)));
end
for k=1:length(S.AB)
    V(k,5) = double(vpa(S.AB(k,1)));
end

% Limit positive real value
n=1;
for k=1:length(V(:,1))
    if isreal(V(k,1)) && isreal(V(k,2)) && isreal(V(k,3)) && isreal(V(k,4)) && isreal(V(k,5))...
            && V(k,1)>0 && V(k,2)>0 && V(k,3)>0 && V(k,4)>0 && V(k,5)>0
        Vreal(n,:) = V(k,:);
        n = n+1;  
    end
end

% Limit heterodimer: Da_TA > && Da_TB >
n=1;
for k=1:length(Vreal(:,1))
    if Da_TA > Vreal(k,5) && Da_TB > Vreal(k,5)
        Vfinal(n,:) = Vreal(k,:);
        n = n+1;  
    end
end


%% Summarize the result
Dtrue_GALL = Vfinal(1);
Dtrue_RALL = Vfinal(2);
Dtrue_GG = Vfinal(3);
Dtrue_RR = Vfinal(4);
Dtrue_GR = Vfinal(5);
Dtrue_G = Dtrue_GALL - 2*Dtrue_GG - Dtrue_GR;
Dtrue_R = Dtrue_RALL - 2*Dtrue_RR - Dtrue_GR;

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