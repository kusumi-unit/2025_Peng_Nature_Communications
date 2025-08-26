clear all

%{
% Solve True T from apparent T
syms Ec Ed Kd Tapp T

% Ec = 3*(Ea+Eb)*(3-Ea-Eb)+6*Ea*Eb
% Ed = (Ea+Eb)^2-2*Ea*Eb
eqn = Tapp == 1/18*(Ec*T-Ed*Kd+Ed*sqrt(Kd^2+6*Kd*T));
S = solve(eqn, T); 
A = S(2);
% T = (3*Ed^2*Kd - Ed*(Kd*(Kd*Ec^2 + 6*Kd*Ec*Ed + 108*Tapp*Ec + 9*Kd*Ed^2))^(1/2) + 18*Ec*Tapp + Ec*Ed*Kd)/Ec^2


% Check the answer
% when Ea =Eb = 1, Ec = 12, Ed = 2
A1 = subs(A, [Ec, Ed], [12, 2]);

% Full equation
syms Ea Eb
Eq1 = 3*(Ea+Eb)*(3-Ea-Eb)+6*Ea*Eb;
Eq2 = (Ea+Eb)^2-2*Ea*Eb;
A2 = subs(A, [Ec, Ed], [Eq1, Eq2]);


% Solve E for [A]
syms E Ea Eb

E = 2*(1-Ea)^2 + 2*(1-Eb)^2 +2*(1-Ea)*(1-Eb) + Ea^2 + Eb^2 +Ea*(1-Eb) + Eb*(1-Ea) + 2*Ea*(1-Ea) + 2*Eb*(1-Eb);
E_ex = expand(E);
%}


% With incidental

syms Aapp1 Bapp1 Aapp2 Bapp2 ABapp Kd EA EB r A Ttrue Tapp

Eq1 = Aapp1 == EA/18*(3*(3-EA-EB)*Ttrue-(EA+EB)*Kd+(EA+EB)*sqrt(Kd^2+6*Kd*Ttrue));
Eq2 = Bapp1 == EB/18*(3*(3-EA-EB)*Ttrue-(EA+EB)*Kd+(EA+EB)*sqrt(Kd^2+6*Kd*Ttrue));
Eq3 = ABapp == EA*EB/18*sqrt(3*Ttrue+Kd-sqrt(Kd^2+6*Kd*Ttrue));

%Eq1 = subs(Eq1, [EA, EB, Kd], [1, 1, 3]);




Eq1_2 = Aapp2 == Aapp1*(2*pi*r^2*(Aapp1+ABapp)/A);
Eq2_2 = Bapp2 == Bapp1*(2*pi*r^2*(Bapp1+ABapp)/A);
Eq1_2 = subs(Eq1_2, lhs(Eq1), rhs(Eq1));
Eq2_2 = subs(Eq2_2, lhs(Eq2), rhs(Eq2));

Eq4 = Tapp == Aapp2 + Bapp2 + 2*ABapp;
Eq4 = subs(Eq4, lhs(Eq1_2), rhs(Eq1_2));
Eq4 = subs(Eq4, lhs(Eq2_2), rhs(Eq2_2));
Eq4 = subs(Eq4, lhs(Eq3), rhs(Eq3));

S = solve(Eq4, Ttrue, 'Real', true);


A1 = subs(S(1), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);
A2 = subs(S(2), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);
A3 = subs(S(3), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);
A4 = subs(S(4), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);
A5 = subs(S(5), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);
A6 = subs(S(6), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);
A7 = subs(S(7), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);
A8 = subs(S(8), [EA, EB, Kd, r, Tapp, A], [1, 1, 3, 0.2, 1, 100]);



