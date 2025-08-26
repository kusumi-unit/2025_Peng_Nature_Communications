function Pos = ErrorAdd2(Pos, Er, PixSize)
for k=1:length(Pos)
    NumSpots = length(Pos{k}(:,1));
    Er_matrix = normrnd(0, Er/PixSize, NumSpots, 2);
    Pos{k}(:,1:2) = Pos{k}(:,1:2) + Er_matrix;
end
end