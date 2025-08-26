clear all; close all;

%% Initial settings
ImSize = 200; % pix x pix square
PixSize = 50; % nm/pix
DtApp = 2; % Apparent total (Red + Green) spots density (copies / um^2)
Nframe = 1;
InciTh = 200; % Threshold distance pf incidental homo dimerization (nm)
FigDraw = 0; % Draw figure of not 
Iteration = 1000;


%% Test
% Initialized the matrix
Nt = DtApp*(ImSize*PixSize/1000)^2;
Nmol_sum = zeros(Nt,Iteration);
for k=1:Iteration
    k
    % Calculate spots #
    Nspots = DtApp*(ImSize*PixSize)^2*10^-6;
    % Generate random tracks
    OriPos = TrackGen(Nframe, Nspots, ImSize);
    % Put # of molecules in the cluster as 0
    OriPos{1}(:,3) = 0;

    [Pos_mod{k}, ~, ~, Nt_ori(k,1), Nt_mod(k,1), Ncluser(k,1), Nmono(k,1)] = InciCalc_forSim(OriPos{1}, PixSize, InciTh, FigDraw);
    
    % Summarize # of molecules
    
    for j=1:length(Pos_mod{k}(:,1))
        Nmol = Pos_mod{k}(j,3);
        Nmol_sum(Nmol,k) = Nmol_sum(Nmol,k)+1;
    end
end

Sum = horzcat(Nt_ori, Nt_mod, Ncluser, Nmono);

Nmol_ave = mean(Nmol_sum,2);

for k=1:length(Nmol_ave(:,1))
    Nmol_proto(k,1) = Nmol_ave(k,1)*k;
end

xxx_Gre_di = sum(Nmol_proto(2:end))/sum(Nmol_proto)*100;
xxx_Gre_tri = sum(Nmol_proto(3:end))/sum(Nmol_proto)*100;
