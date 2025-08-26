function Tracks = TrackGen(Nframe, Nspots, ImSize)
% To make average dimer # in whole movie precisely,
% Copy integer part and decimal part differently
Int_N = floor(Nspots);
Dec_N = Nspots - Int_N;
AddFrame = round(Dec_N*Nframe, 0);
AddFrameIndex = ceil(rand(AddFrame,1)*Nframe);
% Generate Int part
for k=1:Nframe
     Tracks{k}(:,1:2) = rand(Int_N,2) * ImSize;
end
% Add Decimal part
for k=1:AddFrame
    Tracks{AddFrameIndex(k)}(end+1,1:2) = rand(1,2) * ImSize;
end
end