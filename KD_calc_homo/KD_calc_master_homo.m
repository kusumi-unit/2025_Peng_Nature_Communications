clear

% Results for the Gaussian fitting of PCCF 
% y = y0 + (C/(w*sqrt(pi/2)))*exp(-2*((x-xc)/w)^2)
C = 1.55303;
W = 140.67212;

% Experimetal conditions
Area = 100; % observation area, um^2
NG = 39.614; % # of observed green spots / frame
NR = 39.88475; % # of observed red spots / frame
EG = 0.614835511; % Labeling efficiency for green spots: 0-1
ER = 0.663087205; % Labeling efficiency for green spots: 0-1
Th = 200; % threshold distance for incidental colocalizations (nm)


%% Calc
% Convert the unit
A = Area*10^6;

Kd = (((NG + (EG^2/2-2*EG)*( (2*( Th^2/(4*W^2*NG*NR))*( ((A-A*sqrt(1-2*pi*Th^2*NG/A))/(pi*Th^2)) + ((A-A*sqrt(1-2*pi*Th^2*NR/A))/(pi*Th^2)))-( A/(2*pi*NG*NR*W^2))+sqrt(((A/(2*pi*NG*NR*W^2))-2*( Th^2/(4*W^2*NG*NR))*( ((A-A*sqrt(1-2*pi*Th^2*NG/A))/(pi*Th^2)) + ((A-A*sqrt(1-2*pi*Th^2*NR/A))/(pi*Th^2))))^2+12*( Th^2/(4*W^2*NG*NR))*C))/(6*( Th^2/(4*W^2*NG*NR))*EG*ER)) + (pi*Th^2*((A-A*sqrt(1-2*pi*Th^2*NG/A))/(pi*Th^2))^2/(2*A)))/EG)+( (NR + (ER^2/2-2*ER)*( (2*( Th^2/(4*W^2*NG*NR))*( ((A-A*sqrt(1-2*pi*Th^2*NG/A))/(pi*Th^2)) + ((A-A*sqrt(1-2*pi*Th^2*NR/A))/(pi*Th^2)))-( A/(2*pi*NG*NR*W^2))+sqrt(((A/(2*pi*NG*NR*W^2))-2*( Th^2/(4*W^2*NG*NR))*( ((A-A*sqrt(1-2*pi*Th^2*NG/A))/(pi*Th^2)) + ((A-A*sqrt(1-2*pi*Th^2*NR/A))/(pi*Th^2))))^2+12*( Th^2/(4*W^2*NG*NR))*C))/(6*( Th^2/(4*W^2*NG*NR))*EG*ER)) + (pi*Th^2*((A-A*sqrt(1-2*pi*Th^2*NR/A))/(pi*Th^2))^2/(2*A)))/ER))^2/(2*((2*( Th^2/(4*W^2*NG*NR))*( ((A-A*sqrt(1-2*pi*Th^2*NG/A))/(pi*Th^2)) + ((A-A*sqrt(1-2*pi*Th^2*NR/A))/(pi*Th^2)))-( A/(2*pi*NG*NR*W^2))+sqrt(((A/(2*pi*NG*NR*W^2))-2*( Th^2/(4*W^2*NG*NR))*( ((A-A*sqrt(1-2*pi*Th^2*NG/A))/(pi*Th^2)) + ((A-A*sqrt(1-2*pi*Th^2*NR/A))/(pi*Th^2))))^2+12*( Th^2/(4*W^2*NG*NR))*C))/(6*( Th^2/(4*W^2*NG*NR))*EG*ER))*A)*10^6;
