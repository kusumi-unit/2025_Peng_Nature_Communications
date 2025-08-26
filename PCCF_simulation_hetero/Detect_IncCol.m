function Detect_IncCol(Gd,Gm,Rd,Rm,PixSize,Th,DtApp,Kd,Napp_G,Napp_R,Napp_GR,LEr,FilenameRoot)

% for debug
%{
clearvars -except Gd Gm Rd Rm
PixSize = 50; % (nm/pix) 
Th = 200; % Threshold distance (nm)
DtApp = 1;
Kd = 1;
Napp_G = 1;
Napp_R = 1;
Napp_GR = 1;
LEr = 1;
FilenameRoot = 'test';
%}


%% Summarize
Gt = vertcat(Gd,Gm);
Rt = vertcat(Rd,Rm);
% Detect incidentally colocalized spots, and separate
[GmMod,GmInc] = Col(Gm,Rt,Th,PixSize);
[RmMod,RmInc] = Col(Rm,Gt,Th,PixSize);


%% Graph 
close all;
Fig = figure(1);
hold on;
scatter(GmMod(:,1), GmMod(:,2), 20, [0.47,0.67,0.19]);
scatter(RmMod(:,1), RmMod(:,2), 20, [0.85,0.33,0.10]);
scatter(Gd(:,1), Gd(:,2), 40, [0.47,0.67,0.19], 's');
scatter(Rd(:,1), Rd(:,2), 40, [0.85,0.33,0.10], 's');
if isempty(GmInc)
else
    scatter(GmInc(:,1), GmInc(:,2), 40, [0.47,0.67,0.19], '^');
    scatter(RmInc(:,1), RmInc(:,2), 40, [0.85,0.33,0.10], '^');
end
X1 = num2str(DtApp);
X2 = num2str(Kd);
X3 = num2str(Napp_G);
X8 = num2str(Napp_R);
X4 = num2str(Napp_GR);
X6 = num2str(PixSize);
X7 = num2str(LEr);
FigTitle{1} = strcat('Apparent total Density [Red + Green] = ', X1, ' (copies/um^2), Kd = ', X2, ' (copies/um^2)'); 
FigTitle{2} = strcat('Apparent #R:', X8, ', #G:', X3, ', #RG:', X4);
FigTitle{3} = strcat('PixSize = ', X6, ' (nm/pix), GLEr = RLEr = ', X7, ' (nm)');
title(FigTitle);
warning('off', 'MATLAB:legend:IgnoringExtraEntries'); % Legend Ignore alert 
legend('Monomer','', 'Hetero-Dimer', '', 'Incidental homo-Dimer', '')
% Re-scale
xlim([0 200]); ylim([0 200]);
% Fix X and Y axis ratio
pbaspect([1 1 1]);

%%
% Save figure
Filename1 = strcat(FilenameRoot, '.fig');
% Delete the file, if exist
warning('off', 'MATLAB:DELETE:FileNotFound'); % Ignore alert
delete(Filename1);
saveas(Fig, Filename1)
close(1)

% Save to excel file (Xpos, Ypos) in pix unit
Filename2 = strcat(FilenameRoot, '_Spots.xlsx');
% Delete the file, if exist
delete(Filename2);
writecell(num2cell(Gd),Filename2,'Sheet', 'G_di');
writecell(num2cell(GmMod),Filename2,'Sheet', 'G_mono');
writecell(num2cell(GmInc),Filename2,'Sheet', 'G_mono_inc');
writecell(num2cell(Rd),Filename2,'Sheet', 'R_di');
writecell(num2cell(RmMod),Filename2,'Sheet', 'R_mono');
writecell(num2cell(RmInc),Filename2,'Sheet', 'R_mono_inc');


%% subfunction
function [GmMod,GmInc] = Col(Gm,Rt,Th,Psize)
Ninc = 1;
Nmod = 1;
GmInc =[];
for k=1:length(Gm(:,1))
    Dis =[]; % initialize
    for j=1:length(Rt(:,1))
        % Calculate distance
        Dis(j,3) = sqrt((Gm(k,1) - Rt(j,1))^2 + (Gm(k,2) - Rt(j,2))^2);
    end
    if min(Dis) < Th/Psize
        GmInc(Ninc,:) = Gm(k,:);
        Ninc = Ninc+1;
    else
        GmMod(Nmod,:) = Gm(k,:);
        Nmod = Nmod+1;
    end
end


end


end



