function [CalHisto_sum, CalHisto_sum_FL, ColIndex, ColIndex_FL, ...
    Napp_ALL, Napp_GR, Napp_G, Napp_R, ...
    Ntrue_ALL, Ntrue_GR, Ntrue_G, Ntrue_R, ...
    Nreal_GALL, Nreal_RALL, Ninci_G_M, Ninci_R_M,Nreal_GR] = ...
    ColSimu_Gsub2(PixSize, ImSize, DtApp, Kd, LEr, Nframe, BinWidth, ...
    BinMatrix, BinCent, ColIndexCalcN, ColIndexCalcD, NormIndex, ...
    FilenameRoot, Flip, Homo, Inci, InciTh, EG, ER, IndLoop1, IndLoop2)

% Calculate Ture # of hetero-dimers and % of total spots
% Then, reduce the spots by homo-dimerization, labelling efficiency, 
% and incidental colocalization

%% Initial settings for debug

% clear; close all;
% PixSize = 50; % nm/pix
% ImSize = 200; % pix x pix square
% DtApp = 1; % Apparent total (Red + Green) spots density (copies / um^2)
% Kd = 3; % copies / um^2, 0 for stable dimer, inf for random control
% LEr = 100; % Apparent localization error: sigma for single color  
% Nframe = 200; % # of frames for one movie
% BinWidth = 50; % nm
% ColIndexCalcN = [1,2]; % Bin index position for numerator of colocalization index calc
% ColIndexCalcD = [9,10]; % denominator, inf for no normalization
% NormIndex = [inf, inf]; % denominator, and nomalization factor for histogram, inf for no normalization
% EstIndex = [1, 20]; % Bin index position for N dimer estimation 
% FilenameRoot = 'test';
% Flip = 0;
% Homo = 1; % Homo dimer containing or not: 0 or 1
% Inci = 1; % Incidental homo dimer simulation: 1 or 0
% InciTh = 200; % Threshold distance pf incidental homo dimerization (nm)
% EG = 0.7; % Labelling efficiancy for Green: 0~1
% ER = 0.7;
% 
% IndLoop1 = 1;
% IndLoop2 = 1;
% 
% % Bin setting
% BinMax = hypot(ImSize*PixSize, ImSize*PixSize); % Detremine max bin
% BinMatrix = [0:BinWidth:BinMax]; % Create bin matrix
% BinCent(:,1) = BinMatrix(1:end-1) + BinWidth/2; % Bin center matrix


%% Generate random data
% Format for G, Rpos
%[Xpos, Ypos, Dimer ID (1 <= for dimer, 0 for monomer)]

% Calculate the area in um^2 unit
Area = (ImSize * PixSize)^2 * 10^-6;
% Determine # of spots
if Homo == 1
    % Homo-dimer existing case---------------------------------------------
    if Kd ~= inf
        [Dtrue_ALL,Dtrue_GR,Dapp_GR,Dtrue_G,Dtrue_R,Dapp_G,Dapp_R]...
            = SpotCalc_HomoDimer(DtApp,Kd,InciTh,EG,ER,Area);
    else
        [Dtrue_ALL,Dtrue_GR,Dapp_GR,Dtrue_G,Dtrue_R,Dapp_G,Dapp_R]...
            = SpotCalc_NoDimer(DtApp,InciTh,Area);
    end

else
    %Currently not available
    %{
    % No homo-dimer existing case------------------------------------------
    Dtrue_ALL = DtApp;
    Ntrue_ALL = Area*Dtrue_ALL;
    Napp_ALL = Area*DtApp;
    Ntrue_G = Napp_G;
    Ntrue_R = Napp_R;
    Napp_G = ((-Kd+sqrt(Kd^2+2*Kd*DtApp))/2)*Area;
    Napp_R = Napp_G;
    Ntrue_GR = Napp_ALL/2 - NmoApp;
    Napp_GR = Ntrue_GR;
    %}
end

% Convert the density to the N 
Ntrue_ALL = Area*Dtrue_ALL;
Ntrue_GR = Area*Dtrue_GR; 
Napp_GR = Area*Dapp_GR;
Ntrue_G = Area*Dtrue_G;
Ntrue_R = Area*Dtrue_R;
Napp_G = Area*Dapp_G;
Napp_R = Area*Dapp_R;
% Calculate apparent # of spots in a frame
Napp_ALL = 2*Napp_GR + Napp_G + Napp_R;

% Generate randamized tracks
Gpos = TrackGen(Nframe, Ntrue_ALL/2, ImSize);
Rpos = TrackGen(Nframe, Ntrue_ALL/2, ImSize);


%% Post modification 
% Hetero-dimerization 
if Kd == inf
    % No dimer for Kd = inf
    Int_Ndi = 0;
    Int_Ndi_sumG{1} = 0;
    Int_Ndi_sumR{1} = 0;
    % Put all dimer flag as 0
    for k=1:length(Rpos)
        Gpos{k}(:,3) = 0;
        Rpos{k}(:,3) = 0;
    end
else 
    % Copy Gpos data to Rpos data to make dimer
    % To make average dimer # in whole movie precisely,
    % Copy integer part and decimal part differently
    Int_Ndi = floor(Ntrue_GR);
    Dec_Ndi = Ntrue_GR - Int_Ndi;
    AddFrame = round(Dec_Ndi*Nframe, 0);
    % Interger part
    for k=1:length(Rpos)
        Rpos{k}(1:Int_Ndi,:) = Gpos{k}(1:Int_Ndi,:);
        % Add dimer ID: 1 <= for dimer, 0 for monomer
        Rpos{k}(1:Int_Ndi,3) = [1:Int_Ndi];
        Gpos{k}(1:Int_Ndi,3) = [1:Int_Ndi];
    end
    % Decimal part
    for k=1:AddFrame
        Rpos{k}(Int_Ndi+1,:) = Gpos{k}(Int_Ndi+1,:);
        Rpos{k}(Int_Ndi+1,3) = Int_Ndi+1;
        Gpos{k}(Int_Ndi+1,3) = Int_Ndi+1;
    end
end

% Add position determination error
Gpos = ErrorAdd2(Gpos, LEr, PixSize);
Rpos = ErrorAdd2(Rpos, LEr, PixSize);

% Homo-dimerization 
Gpos = HomoGen(Gpos,Ntrue_GR/2,Nframe);
Rpos = HomoGen(Rpos,Ntrue_GR/2,Nframe);

% Delete unlabelled spots rondomly
Gpos = DelUnlab(Gpos,EG,Nframe);
Rpos = DelUnlab(Rpos,ER,Nframe);

% Delete incidental dimer
if Inci == 1
    for k=1:Nframe
        % Skip if < 1 spots
        if length(Gpos{k}(:,1)) > 1
            [Gpos_mod{k},~,~,Nt_oriG(k),Nt_modG(k),~,~] = InciCalc2(Gpos{k}, PixSize, InciTh, 0);
            Ninci_G(k) = Nt_oriG(k) - Nt_modG(k);
        else
            Gpos_mod{k} = Gpos{k};
        end
        if length(Rpos{k}(:,1)) > 1
            [Rpos_mod{k},~,~,Nt_oriR(k),Nt_modR(k),~,~] = InciCalc2(Rpos{k}, PixSize, InciTh, 0);
            Ninci_R(k) = Nt_oriR(k) - Nt_modR(k);
        else
            Rpos_mod{k} = Rpos{k};
        end
    end
    Ninci_G_M = mean(Ninci_G);
    Ninci_R_M = mean(Ninci_R);  
else
    Gpos_mod = Gpos;
    Rpos_mod = Rpos;
    Ninci_G_M = 0; Ninci_R_M = 0;
end

% Delete un-paired dimer ID
Gpos_mod = CorrDimerID(Gpos_mod,Rpos_mod,Nframe);
Rpos_mod = CorrDimerID(Rpos_mod,Gpos_mod,Nframe);

% Sort by dimer ID
for k=1:Nframe
    [~,I] = sort(Gpos_mod{k}(:,3),'descend');
    Gpos_mod{k} = Gpos_mod{k}(I,:);
end
for k=1:Nframe
    [~,I] = sort(Rpos_mod{k}(:,3),'descend');
    Rpos_mod{k} = Rpos_mod{k}(I,:);
end

% Count % of spots
[Nreal_GALL,Nreal_GR,Ndi_1st] = SpotsStat(Gpos_mod,Nframe);
[Nreal_RALL,~,~] = SpotsStat(Rpos_mod,Nframe);
Dreal_ALL = (Nreal_GALL+Nreal_RALL)/Area;


%% Graph
% Generate Fig for only the first loop
if IndLoop1 == 1 && IndLoop2 == 1
    % Summarize data
    Gm(:,1) = Gpos_mod{1}(Ndi_1st+1:end,1);
    Gm(:,2) = Gpos_mod{1}(Ndi_1st+1:end,2);
    Rm(:,1) = Rpos_mod{1}(Ndi_1st+1:end,1);
    Rm(:,2) = Rpos_mod{1}(Ndi_1st+1:end,2);
    Gd(:,1) = Gpos_mod{1}(1:Ndi_1st,1);
    Gd(:,2) = Gpos_mod{1}(1:Ndi_1st,2);
    Rd(:,1) = Rpos_mod{1}(1:Ndi_1st,1);
    Rd(:,2) = Rpos_mod{1}(1:Ndi_1st,2);    
    % for incidental col detection for figure  
    Detect_IncCol(Gd,Gm,Rd,Rm,PixSize,InciTh,DtApp,Kd,Napp_G,Napp_R,Napp_GR,LEr,FilenameRoot)
end

% Just in case, delete excess ID data
for k=1:Nframe
    Gpos_mod{k}(:,3) = [];
    Rpos_mod{k}(:,3) = [];
end


%% Distance calc
%  Create Calibration curve

% % Line calibration
% Sr = Lr_func(BinCent, Nreal_GALL, Nreal_RALL, Nframe, BinWidth, Area);

% Naito's Sr curve
[Sr] = Sr_func(BinWidth, ImSize, ImSize, PixSize, BinMatrix, Nreal_GALL/Area, Nreal_RALL/Area, Nframe);
% Because the calibration curve start from 0, delete head value to match count histogram
Sr = Sr(2:end, 1);
% Data summary
BinCent_Norm = BinCent;

%  Create Raw data histogram
Histo_Raw = COL2(Gpos_mod, Rpos_mod, PixSize, BinMatrix); 
[CalHisto_sum, ColIndex] = StatSum(Histo_Raw, Sr, NormIndex, ColIndexCalcN, ColIndexCalcD);

%  Create flipped data and histogram
if Flip == 1
    for k=1:length(Gpos_mod)
        dx = (Gpos_mod{k}(:,1) - ImSize/2); % Calculate distance to center
        dy = (Gpos_mod{k}(:,2) - ImSize/2);
        Gpos_Fl{k}(:,1) = Gpos_mod{k}(:,1) - 2*dx; % Parallele shift to the center
        Gpos_Fl{k}(:,2) = Gpos_mod{k}(:,2) - 2*dy;
    end
    % Calculate histogram
    Histo_FL = COL2(Gpos_Fl, Rpos_mod, PixSize, BinMatrix);
    [CalHisto_sum_FL, ColIndex_FL] = StatSum(Histo_FL, Sr, NormIndex, ColIndexCalcN, ColIndexCalcD);
else
    Histo_FL = [];
    CalHisto_sum_FL = 0; ColIndex_FL = 0;
end


%% Sub functions
function [CalHisto_sum, ColIndex] = StatSum(Histo, Sr, NormIndex, ColIndexCalcN, ColIndexCalcD) 
% Normalized by Naito function
CalHisto_sum = Histo ./ Sr;
% Re-nomalize by Index position
if NormIndex(1) == inf && NormIndex(2) == inf
    % Skip normalization for inf input
else
    HistoNormFactor = mean(CalHisto_sum(NormIndex(1):NormIndex(2),1));
    CalHisto_sum = CalHisto_sum/HistoNormFactor;
end
% Calculate colocalization index
if ColIndexCalcD(1) == inf && ColIndexCalcD(2) == inf
    % Skip normalization
    ColIndex = mean(CalHisto_sum(ColIndexCalcN(1):ColIndexCalcN(2),1));
else
    ColNormFactor = mean(CalHisto_sum(ColIndexCalcD(1):ColIndexCalcD(2),1));
    ColIndex = mean(CalHisto_sum(ColIndexCalcN(1):ColIndexCalcN(2),1))/ColNormFactor;
end

end

function Pos = HomoGen(Pos,Ntrue_Homo,Nframe)
Int = floor(Ntrue_Homo);
Dec = Ntrue_Homo - Int;
AddF = round(Dec*Nframe, 0);
% Interger + Decimal part
for H=1:Nframe
    % if the frame# is 1 ~ Addfrane, # of homodimer is Int_Ndi+1 
    if H <= AddF
        L = Int+1;
    else
        L = Int;
    end
    E = length(Pos{H}(:,1)); % determine the end
    % Duplicate data from the end
    Pos{H}(E-L+1:E,:) = Pos{H}(E-2*L+1:E-L,:);
    % Put homo-dimer flag 
    Pos{H}(E-2*L+1:E,3) = 0; % Currently not using for homodimer
end
end

function Pos = DelUnlab(Pos,E,Nframe)
% Count total # of spots
for D=1:Nframe
    Lsum(D) = length(Pos{D}(:,1));
end
Ntotal = sum(Lsum);
% Calculate # to delete
Ndel = round(Ntotal*(1-E));
% Delete randomly selected spots
for D=1:Ndel
    % select frame and pos index
    % Re-select the frame if empty
    DelFrame = randi(Nframe,1);
    while isempty(Pos{DelFrame})
       DelFrame = randi(Nframe,1);
    end
    DelIdx = randi(length(Pos{DelFrame}(:,1)),1);
    % delete
    Pos{DelFrame}(DelIdx,:) = [];
end
end

function Pos1 = CorrDimerID(Pos1,Pos2,Nframe)
for k1=1:Nframe
    if isempty(Pos1{k1}) || isempty(Pos2{k1}) % skip empty frame
    else
        for j1=1:length(Pos1{k1}(:,1))
            if Pos1{k1}(j1,3) == 0
            else
                % find paired hetero-dimer
                Flag1 = [];
                Flag1 = Pos1{k1}(j1,3) == Pos2{k1}(:,3);
                % If no-pair find, put ID as 0
                if sum(Flag1) == 0
                    Pos1{k1}(j1,3) = 0;  
                end
            end
        end
    end
end
end

function [Nt,Ndi,Ndi_1st] = SpotsStat(Pos,Nframe)
for k2=1:Nframe
    % Calculate total # of spots 
    Tot(k2) = length(Pos{k2}(:,1));
    % Calculate total # of dimer
    Flag2 = Pos{k2}(:,3) == 0;
    NDim(k2) = Tot(k2) - sum(Flag2);
end
Nt = sum(Tot)/Nframe;  
Ndi = sum(NDim)/Nframe;
Ndi_1st = NDim(1);
end


end