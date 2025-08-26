function Raw_Histo = COL2(Gpos, Rpos, PixSize, BinMatrix)
% Calculate # id total distance
for k=1:length(Gpos)
    Dis(k) = length(Gpos{k}(:,1)) * length(Rpos{k}(:,1));
end
NumTotal = sum(Dis);
% Generate inititial sum matrix
Dis_sum = nan(NumTotal,1);
n = 1;
for k=1:length(Gpos) % Loop for frame#    
for i=1:length(Gpos{k}(:,1)) % Loop for Gpos
    % Calaculate X distance between Gpos(i) and Rpos(:) and input in Colomn3 of Rpos
    % Y distance to Colomn4, R distance to Colomn5
    Rpos{k}(:,3) = Gpos{k}(i,1) - Rpos{k}(:,1);         
    Rpos{k}(:,4) = Gpos{k}(i,2) - Rpos{k}(:,2);         
    Rpos{k}(:,5) =  hypot(Rpos{k}(:,3), Rpos{k}(:,4));
    L = length(Rpos{k}(:,5));
    Dis_sum(n:n+L-1, 1) = Rpos{k}(:,5); % Summarize data
    n = n+L;
end
end
% Re-scale
Dis_sum = Dis_sum*PixSize;
% Create histogram
figure(2) % Fig(1) for GUI
Hist = histogram(Dis_sum, BinMatrix); 
Raw_Histo(:,1) = Hist.BinCounts;
Norm_Histo(:,1) = Raw_Histo/sum(Raw_Histo);
close(2)
end