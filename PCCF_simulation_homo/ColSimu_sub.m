function ColSimu_sub(PixSize, ImSize, Dapp_ALL, Kd, LEr, Nframe, ...
BinWidth, ColIndexCalcN, ColIndexCalcD, NormIndex, GaussFitIndex, FixYoff, ...
Flip, Homo, Inci, InciTh, EG, ER, NumSample, Iteration, Loop)

% Calculate Ture # of hetero-dimers and % of total spots
% Then, reduce the spots by homo-dimerization, labelling efficiency, 
% and incidental colocalization

% curve fitting toolbox required

%% Initial settings for debug
% clear; close all;
% PixSize = 50; % nm/pix
% ImSize = 200; % pix x pix square
% Dapp_ALL = 1; % Apparent total (Red + Green) spots density (copies / um^2)
% Kd = 3; % copies / um^2, 0 for stable dimer, inf for random control
% LEr = 100; % Apparent localization error: sigma for single color  
% Nframe = 200; % # of frames for one movie
% BinWidth = 50; % nm
% ColIndexCalcN = [1,2]; % Bin index position for numerator of colocalization index calc
% ColIndexCalcD = [9,10]; % denominator, inf for no normalization
% NormIndex = [inf,inf]; % denominator, and nomalization factor for histogram, inf for no normalization
% FilenameRoot = 'test';
% GaussFitIndex = 20; % Bin index position for Gauss fitting range 1~
% FixYoff = 0; % Fix Y0 as 1 or not for Gauss fitting  
% Flip = 1;
% Homo = 1; % Homo dimer simulation: 1 or 0
% Inci = 1; % Incidental homo dimer simulation: 1 or 0
% InciTh = 200; % Threshold distance pf incidental homo dimerization (nm)
% EG = 0.7; % Labelling efficiancy for Green: 0~1   0.74, 0.7
% ER = 0.7;
% NumSample = 20;
% Iteration = 1;
% 
% Loop = 1;


%% Calculation
% Calculate the area in um^2 unit
Area = (ImSize * PixSize)^2 * 10^-6;
% Bin setting
BinMax = hypot(ImSize*PixSize, ImSize*PixSize); % Detremine max bin
BinMatrix = [0:BinWidth:BinMax]; % Create bin matrix
BinCent(:,1) = BinMatrix(1:end-1) + BinWidth/2; % Bin center matrix

% Create result file name
FilenameRoot = strcat('AppDens-', num2str(Dapp_ALL), '_Kd-', num2str(Kd), ...
'_Homo-', num2str(Homo), '_Inci-', num2str(Inci), '_L', num2str(Loop));
for m=1:Iteration
    CurrentLoop_CurrentIter = [Loop, m]
    % Initialization
    Histo=[]; Histo_FL =[]; ColInd=[]; ColInd_FL=[];
    Napp_ALL=[]; Napp_GR=[]; Napp_G=[]; Napp_R=[]; 
    Ntrue_ALL=[]; Ntrue_GR=[]; Ntrue_G=[]; Ntrue_R=[];
    Nreal_GALL=[]; Nreal_RALL=[]; Ninci_G=[]; Ninci_R=[]; Nreal_GR=[];
    parfor k=1:NumSample
        [Histo(:,k), Histo_FL(:,k), ColInd(k), ColInd_FL(k), ...
        Napp_ALL(k), Napp_GR(k), Napp_G(k), Napp_R(k), ...
        Ntrue_ALL(k), Ntrue_GR(k), Ntrue_G(k), Ntrue_R(k), ...
        Nreal_GALL(k), Nreal_RALL(k), Ninci_G(k), Ninci_R(k),Nreal_GR(k)] = ...
         ColSimu_Gsub2(PixSize, ImSize, Dapp_ALL, Kd, LEr, Nframe, BinWidth, ...
         BinMatrix, BinCent, ColIndexCalcN, ColIndexCalcD, NormIndex, FilenameRoot, ...
         Flip, Homo, Inci, InciTh, EG, ER, m, k);
    end   
    [ColInd_M(m,1), ColInd_SE(m,1), Histo_M(:,m), Histo_SE(:,m)] = STAT1(ColInd, Histo);
    
    if Flip == 1
        [ColInd_M_FL(m,1), ColInd_SE_FL(m,1), Histo_M_FL(:,m), Histo_SE_FL(:,m)] = STAT1(ColInd_FL, Histo_FL);
    else
        ColInd_M_FL = []; ColInd_SE_FL = []; Histo_M_FL = []; Histo_SE_FL = [];
    end
    [EstNab(m,:),EstNa(m,:),EstNb(m,:),FitResults(m,:),EstKd(m,1),EstKd_SE(m,1),...
    NinciG_M(m,1),NinciG_SE(m,1),NinciR_M(m,1),NinciR_SE(m,1),NtRealG_M(m,1),NtRealR_M(m,1),Nreal_GR_M(m,1)] = ...
    STAT2(Napp_GR(1),Napp_G(1),Napp_R(1),Nreal_GALL,Nreal_RALL,Area,LEr,Histo_M(:,m),BinCent,...
    GaussFitIndex,Homo,EG,ER,Ninci_G,Ninci_G,FixYoff,Inci,InciTh,Nreal_GR);
end

%(NtRealG_M(1) + NtRealR_M(1))/Area


%% Data summary
% Determine the file name
Filename = strcat(pwd,'/Results/', FilenameRoot, '.xlsx');
% Delete the file, if exist
warning('off', 'MATLAB:DELETE:FileNotFound') % ignore alert
delete(Filename);
% Prepare Basic information
Sum1{1,1} = 'PixSize (nm)';
Sum1{2,1} = 'ImSize (pix x pix square)';
Sum1{3,1} = 'Apparent total density (spots / um^2) for seeding';
Sum1{4,1} = 'True total density (copies / frame)';
Sum1{5,1} = 'Real total G density (spots / frame)';
Sum1{6,1} = 'Real total R density (spots / frame)';
Sum1{7,1} = 'Kd (copies / um^2)';
Sum1{8,1} = 'Labeling efficiency for G';
Sum1{9,1} = 'Labeling efficiency for R';
Sum1{10,1} = 'True # of G monomer / frame';
Sum1{10,1} = 'True # of G monomer / frame';
Sum1{11,1} = 'Mean estimated true # of G monomer / frame';
Sum1{12,1} = 'Apparent # of G monomer / frame';
Sum1{13,1} = 'True # of R monomer / frame';
Sum1{14,1} = 'Mean estimated true # of R monomer / frame';
Sum1{15,1} = 'Apparent # of R monomer / frame';
Sum1{16,1} = 'True # of GR dimer / frame';
Sum1{17,1} = 'Mean estimated true # of GR dimer / frame';
Sum1{18,1} = 'Apparent # of GR dimer / frame';
Sum1{19,1} = 'Real # of GR dimer / frame';
Sum1{20,1} = 'Localization error: sigma (nm) *red, green each ';
Sum1{21,1} = 'Range for numerator of colocalization index (nm)';
Sum1{22,1} = 'Range for denominator of colocalization index (nm)';
Sum1{23,1} = 'Range for the histogram normalization (nm)';
Sum1{24,1} = 'Range for the Gaussian fitting (nm)';
Sum1{25,1} = 'Fix Gauss offset';
Sum1{26,1} = '# of frames (per movie)';
Sum1{27,1} = '# of movies (per dataset)';
Sum1{28,1} = 'Iteration';
Sum1{29,1} = 'Homo-dimer simulation';
Sum1{30,1} = 'Incidental homo-dimer simulation';
Sum1{31,1} = 'Incidental homo-dimer threshold (nm)';
Sum1{1,2} = num2str(PixSize);
Sum1{2,2} = num2str(ImSize);
Sum1{3,2} = num2str(Dapp_ALL);
Sum1{4,2} = num2str(Ntrue_ALL(1));
Sum1{5,2} = num2str(mean(NtRealG_M));
Sum1{6,2} = num2str(mean(NtRealR_M));
Sum1{7,2} = num2str(Kd);
Sum1{8,2} = num2str(EG);
Sum1{9,2} = num2str(ER);
Sum1{10,2} = num2str(Ntrue_G(1));
Sum1{11,2} = num2str(mean(EstNa));
Sum1{12,2} = num2str(Napp_G(1));
Sum1{13,2} = num2str(Ntrue_R(1));
Sum1{14,2} = num2str(mean(EstNb));
Sum1{15,2} = num2str(Napp_R(1));
Sum1{16,2} = num2str(Ntrue_GR(1));
Sum1{17,2} = num2str(mean(EstNab));
Sum1{18,2} = num2str(Napp_GR(1));
Sum1{19,2} = num2str(mean(Nreal_GR_M));
Sum1{20,2} = num2str(LEr);
Sum1{21,2} = IgnInf(BinCent, BinWidth, ColIndexCalcN);
Sum1{22,2} = IgnInf(BinCent, BinWidth, ColIndexCalcD);
Sum1{23,2} = IgnInf(BinCent, BinWidth, NormIndex);
Sum1{24,2} = strcat('0 - ', num2str(BinCent(GaussFitIndex)));
if FixYoff==1
    Sum1{25,2} = 'Yes';
else
    Sum1{25,2} = 'No';
end
Sum1{26,2} = num2str(Nframe);
Sum1{27,2} = num2str(NumSample);
Sum1{28,2} = num2str(Iteration);
if Homo==1
    Sum1{29,2} = 'Yes';
else
    Sum1{29,2} = 'No';
end
if Inci==1
    Sum1{30,2} = 'Yes';
else
    Sum1{30,2} = 'No';
end
Sum1{31,2} = num2str(InciTh);

% Data summary
Sum2 = cell(Iteration+5,15);
Sum2{1,1} = 'Raw'; Sum2{1,7} = 'y = C * exp(-x^2/2w^2) + y0'; Sum2{1,16} = 'Flipped';  Sum2{1,19} = '# of incidental homo-dimer (in a frame)'; 
Sum2(2,1:22) = {'Col Index', '', '', 'Estimated Kd', '', '', 'C', '', '', 'w', '', '', 'y0', '', '', 'Col Index', '', '', 'G', '', '', 'R'};
Sum2(3,1:23) = {'Mean of Mean', 'Mean of SE', '', 'Mean of Mean', 'Mean of SE', '', 'Mean of Mean', 'Mean of SE', '', 'Mean of Mean', 'Mean of SE', '', 'Mean of Mean', 'Mean of SE', '', 'Mean of Mean', 'Mean of SE', '', 'Mean of Mean', 'Mean of SE', '', 'Mean of Mean', 'Mean of SE'};
Sum2(5,1:23) = {'Mean', 'SE', '', 'Mean', 'SE', '', 'Value', 'SE', '', 'Value', 'SE', '', 'Value', 'SE', '', 'Mean', 'SE', '', 'Mean', 'SE', '', 'Mean', 'SE'};

Sum2{4,1} = num2str(mean(ColInd_M));
Sum2{4,2} = num2str(mean(ColInd_SE));
Sum2{4,4} = num2str(mean(EstKd));
Sum2{4,5} = num2str(mean(EstKd_SE));
Sum2{4,7} = num2str(mean(FitResults(:,1)));
Sum2{4,8} = num2str(mean(FitResults(:,2)));
Sum2{4,10} = num2str(mean(FitResults(:,3)));
Sum2{4,11} = num2str(mean(FitResults(:,4)));
Sum2{4,13} = num2str(mean(FitResults(:,5)));
Sum2{4,14} = num2str(mean(FitResults(:,6)));
Sum2{4,16} = num2str(mean(ColInd_M_FL));
Sum2{4,17} = num2str(mean(ColInd_SE_FL));
Sum2{4,19} = num2str(mean(NinciG_M));
Sum2{4,20} = num2str(mean(NinciG_M));
Sum2{4,22} = num2str(mean(NinciR_M));
Sum2{4,23} = num2str(mean(NinciR_M));
Sum2(6:5+Iteration,1) = num2cell(ColInd_M);
Sum2(6:5+Iteration,2) = num2cell(ColInd_SE);
Sum2(6:5+Iteration,4) = num2cell(EstKd);
Sum2(6:5+Iteration,5) = num2cell(EstKd_SE);
Sum2(6:5+Iteration,7) = num2cell(FitResults(:,1));
Sum2(6:5+Iteration,8) = num2cell(FitResults(:,2));
Sum2(6:5+Iteration,10) = num2cell(FitResults(:,3));
Sum2(6:5+Iteration,11) = num2cell(FitResults(:,4));
Sum2(6:5+Iteration,13) = num2cell(FitResults(:,5));
Sum2(6:5+Iteration,14) = num2cell(FitResults(:,6));
if Flip == 1
Sum2(6:5+Iteration,16) = num2cell(ColInd_M_FL);
Sum2(6:5+Iteration,17) = num2cell(ColInd_SE_FL);
end
Sum2(6:5+Iteration,19) = num2cell(NinciG_M);
Sum2(6:5+Iteration,20) = num2cell(NinciG_M);
Sum2(6:5+Iteration,22) = num2cell(NinciR_M);
Sum2(6:5+Iteration,23) = num2cell(NinciR_M);

% Save to xls file
warning('off', 'MATLAB:xlswrite:AddSheet'); % Ignore alert 
writecell(Sum1,Filename,'Sheet', 'Parameters');
writecell(Sum2,Filename,'Sheet', 'Summary');
% Histogram output
DataSum(Iteration, BinCent, Histo_M, Filename, 'Histo_Mean')
DataSum(Iteration, BinCent, Histo_SE, Filename, 'Histo_SE')
if Flip == 1
DataSum(Iteration, BinCent, Histo_M_FL, Filename, 'FL Histo_Mean')
DataSum(Iteration, BinCent, Histo_SE_FL, Filename, 'FL Histo_SE')
end
% Delete extra sheet
%DeleteXlsSheet(Filename, [1 3]);


%% sub functions
function [ColInd_M, ColInd_SE, Histo_M, Histo_SE] = STAT1(ColInd, Histo)
ColInd_M = mean(ColInd);
ColInd_SD = std(ColInd);
ColInd_SE = ColInd_SD / sqrt(length(ColInd));
Histo_M = mean(Histo, 2);
Histo_SD = nanstd(Histo,0,2);
Histo_SE = Histo_SD / sqrt(length(ColInd));
end

function [EstNab,EstNa,EstNb,FitResults,EstKd,EstKd_SE,NinciG_M,NinciG_SE,...
          NinciR_M,NinciR_SE,NtRealG_M,NtRealR_M,Nreal_GR_M] = ...
          STAT2(Napp_GR,Napp_G,Napp_R,Nreal_GALL,Nreal_RALL,Area,LEr,Histo_M,BinCent,...
          GaussFitIndex,Homo,EG,ER,Ninci_G,Ninci_R,FixYoff,Inci,InciTh,Nreal_GR)
% Calculate % deleted as incidental homo-dimer
NinciG_M = mean(Ninci_G);
NinciG_SE = std(Ninci_G);
NinciR_M = mean(Ninci_R);
NinciR_SE = std(Ninci_R);
Nreal_GR_M = mean(Nreal_GR);
% Calculate # of counted spots 
NtRealG_M = mean(Nreal_GALL);
NtRealR_M = mean(Nreal_RALL);
% Gaussian fitting
% Determine inicial parameter
P = (Napp_G+Napp_GR)*(Napp_R+Napp_GR);
IniY0 = (P-Napp_GR)/P;
IniW = sqrt(2*(LEr^2));
IniC = Area*10^6*Napp_GR/(P*2*pi*IniW^2); % convert Area unit um to nm,
% Fitting: [C * exp(-(x^2)/(2*W^2)) + 1] or [C * exp(-(x^2)/(2*W^2)) + Y0]
if FixYoff == 1
    GaussResults = GaussFit3_2Para(BinCent(1:GaussFitIndex), Histo_M(1:GaussFitIndex), IniC, IniW);
    GaussResults(3,1) = 1; GaussResults(3,2) = 0; 
else
    GaussResults = GaussFit3(BinCent(1:GaussFitIndex), Histo_M(1:GaussFitIndex), IniC, IniW, IniY0);
end
% Output format: C, C_SE; W, W_SE; Y0, Y0_SE;
% Calculate parameters
[EstNab,EstNa,EstNb,EstKd,EstKd_SE,FitResults] = ...
    ParaEst(Area,NtRealG_M,NtRealR_M,EG,ER,InciTh,GaussResults,Inci,Homo);

end


function Out = IgnInf(BinCent, BinWidth, Input)
if Input(1) == inf || Input(2) == inf
    % Ignore inf input
    Out = 'NA';
else
    Out = strcat(num2str(BinCent(Input(1))-BinWidth/2), ' - ', num2str(BinCent(Input(2))+BinWidth/2));
end
end

function DataSum(Iteration, BinCent, Histo, Filename, SheetName)
Header{1,1} = 'Bin center (nm)';
Header{1,Iteration+3} = 'Mean';
BinCent_sum = num2cell(BinCent(:,1));
Histo_M_sum = num2cell(Histo);
Blank = cell(length(BinCent(:,1)),1);
Mean = num2cell(mean(Histo, 2));
Sum = horzcat(BinCent_sum, Histo_M_sum, Blank, Mean);
Sum = vertcat(Header, Sum);
writecell(Sum, Filename, 'Sheet', SheetName);
end

function DeleteXlsSheet(Filename, SheetIndex)
% The xls file should be in the current folder
% XLSINFO を使用してワークブックの情報を取得
XL_file = [pwd strcat('\', Filename)];
[type, sheet_names] = xlsfinfo(XL_file);
% Excel を COM オートメーションサーバとしてオープン
Excel = actxserver('Excel.Application');
% アプリケーションを隠す
set(Excel, 'Visible', 0);
% エクセルのアラートを非表示に設定
set(Excel,'DisplayAlerts',0);
% エクセルのワークブックのハンドルを取得
Workbooks = Excel.Workbooks;
% エクセルワークブックを開き、アクティベート
Workbook = Workbooks.Open(XL_file);
% アクティブなワークブックの中からシートを選択
Sheets = Excel.ActiveWorkBook.Sheets;
index_adjust = 0;
% ユーザの入力に合わせてワークシートを削除
for i = SheetIndex(1):SheetIndex(2)
      current_sheet = get(Sheets, 'Item', (i-index_adjust));
      invoke(current_sheet, 'Delete')
      out_string = sprintf('ワークシート %s は削除されました。',sheet_names{i});
      index_adjust = index_adjust +1;
end
% ワークブックを保存
Workbook.Save;
% ワークブックをクローズ
Workbooks.Close;
% エクセルを終了
invoke(Excel, 'Quit');
% ActiveX オブジェクトのハンドルを消去
delete(Excel);
end


end