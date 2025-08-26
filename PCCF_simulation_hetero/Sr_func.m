function [NCI] = Sr_func(Bin, XSize, YSize, Scale, BinMatrix, DensityG, DensityR, NumFrame)
%{
Function for Col_Index
Calculate the calibration curve for pair-wise distance and histogram
Based on the Naito-kun's theory:
f(r) = interval * (r - interval / 2) * (2* 4 * pi() * x * scale * y * ...
scale - 8 * r * (x * scale + y * scale) / 2 + 2 * r ^ 2)

Note dat NCI(R) means area density of "R to R - interval" dounut area.

180313 Taka Tsunoyama
%}
Naito(1) =0;
for k=2:length(BinMatrix)
    R = BinMatrix(k);
    Naito(k) = Bin*(R-Bin/2)*(2*pi*XSize*Scale*YSize*Scale-8*R*(XSize*Scale+YSize*Scale)/2+2*R^2);
    % if value is < 0, put 0 
    if Naito(k) < 0
        Naito(k) =0;
    end
end
% Normalized by density
% Total density in a frame: convert /um^2 to /nm^2
DensityG = DensityG * (10^-6);
DensityR = DensityR * (10^-6);
NormFactor = (DensityG*DensityR) * NumFrame; % Multiplied by frame#
NCI(:,1) = Naito * NormFactor; % convert % to ratio

%NCI(:,1) = Naito/sum(Naito);
end