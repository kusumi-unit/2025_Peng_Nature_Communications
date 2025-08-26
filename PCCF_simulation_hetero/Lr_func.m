function Lr = Lr_func(BinCent, NG, NR, Nframe, BinWidth, Area)

Lr = 2*pi*BinCent;
Area_nm = Area*10^6;
NormFactor = BinWidth*NG*NR*Nframe/Area_nm;
Lr = Lr * NormFactor;

end