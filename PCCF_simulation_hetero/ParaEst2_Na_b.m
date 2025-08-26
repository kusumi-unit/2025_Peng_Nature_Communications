function [EstNa,EstNb,EstNaa,EstNbb,Estkon3,EstKd]=...
    ParaEst2_Na_b(EstNab,NtRealG,NtRealR,DobsG,DobsR,...
    kon1,koff1,kon2,koff2,koff3,InciTh,EG,ER,A)

% G + G <--> GG or A + A <--> AA: kon1, koff1
% R + R <--> RR or B + B <--> BB: kon2, koff2
% G + R <--> GR or A + B <--> AB: kon3, koff3

%{
% for debug
clear all;
% True Kd: 3
% True kon3: 2.666
% True GR: 11.4616
% True kon3: 
% Ntrue G all: 75.30
% Ntrue R all: 75.3
% Ntrue GG: 2.11;
% Ntrue RR: 4.86;
%{
kon1 = 0.47; % homo-DOR before stimlation
koff1 = 8.00;
kon2 = 1.03; % homo-KOR before stimulation
koff2 = 6.71;
koff3 = 8; % hetero-dimer
%}
Comb = 'MD'; % MD or KD
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

A = 100;
NtRealG = 0.05;
NtRealR = 0.05;
DobsG = 0.05;
DobsR = 0.05;
EstNab = 50; %11.4616;  % true value 
EG = 0.7;
ER = 0.7;
InciTh = 200;
Inci = 1;
Homo = 1;
%}


%% Calculate Na and Nb

%--------------------------------------------------------------------------
% for simulation, use the value simulated seed parameter gives better results
% Disable this part for the calculation with real density
% NtRealG = DobsG*A;
% NtRealR = DobsR*A;
%--------------------------------------------------------------------------

% Convert the unit
Th = InciTh/1000;
EA = EG; EB = ER;

EstDAB = EstNab/A;
Napp_Gall = (A-A*sqrt(1-2*pi*Th^2*NtRealG/A))/(pi*Th^2);
Napp_Rall = (A-A*sqrt(1-2*pi*Th^2*NtRealR/A))/(pi*Th^2);
Ninc_TG = pi*Th^2*Napp_Gall^2/(2*A);
Ninc_TR = pi*Th^2*Napp_Rall^2/(2*A);

Dobs_TA = NtRealG/A;
Dobs_TB = NtRealR/A;
% TA: true density of toal A 
syms TA TB AA BB kon3
eqn1 = kon1*(TA-2*AA-EstDAB)^2-koff1*AA == 0;
eqn2 = kon2*(TB-2*BB-EstDAB)^2-koff2*BB == 0;
eqn3 = kon3*(TA-2*AA-EstDAB)*(TB-2*BB-EstDAB)-koff3*EstDAB == 0;
eqn4 = Dobs_TA == EA*(TA-2*AA-EstDAB) + EA*(2-EA)*AA + EA*EstDAB - Ninc_TG/A;
eqn5 = Dobs_TB == EB*(TB-2*BB-EstDAB) + EB*(2-EB)*BB + EB*EstDAB - Ninc_TR/A;
% Solve functions
S = solve([eqn1,eqn2,eqn3,eqn4,eqn5], [TA,TB,AA,BB,kon3]);

% Summarize the result
% V = [TA, TB, AA, BB, kon3]
V= [];
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
for k=1:length(S.kon3)
    V(k,5) = double(vpa(S.kon3(k,1)));
end

% Limit positive real value
n=1;
Vreal= [];
if isempty(V)
else
    for k=1:length(V(:,1))
        if isreal(V(k,1)) && isreal(V(k,2)) && isreal(V(k,3)) && isreal(V(k,4)) && isreal(V(k,5))...
                && V(k,1)>0 && V(k,2)>0 && V(k,3)>0 && V(k,4)>0 && V(k,5)>0
            Vreal(n,:) = V(k,:);
            n = n+1;  
        end
    end
end
% Limit # of dimer: Dobs_TA/2 > AA && Da_TB/2 > BB && Da_TA+Da_TB > AB*2
n=1;
Vfinal= [];
if isempty(Vreal)
    else
    for k=1:length(Vreal(:,1))
        if Vreal(k,1) > Vreal(k,3)*2+EstDAB &&...
           Vreal(k,2) > Vreal(k,4)*2+EstDAB &&...
           Vreal(k,1)+Vreal(k,2) > EstDAB*2
            Vfinal(n,:) = Vreal(k,:);
            n = n+1;  
        end
    end
end

% Avoid null values
if isempty(Vfinal)
    Vfinal = [nan,nan,nan,nan,nan];
end


%% Summarize the result
EstNa = (Vfinal(1)-2*Vfinal(3)-EstNab/A)*A;
EstNb = (Vfinal(2)-2*Vfinal(4)-EstNab/A)*A;
EstNaa = Vfinal(3)*A;
EstNbb = Vfinal(4)*A;
Estkon3 = Vfinal(5);
EstKd = koff3/Estkon3;


end