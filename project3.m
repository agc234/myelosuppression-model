%% Part 1: System Analysis of the Friberg's Model
% The goal is to plot neutrophils vs days under Edrugss = 0.250
% and to find the Edrugss what ANC = 1.5.
circ0 = 5.04;
ANC = 1.5;
gamma = 0.16;
kcirc = 4/120.4;
ktr = 4/120.4;
kprol = 4/120.4;
alpha = 1;
alpha_max = circ0/ANC;
Edrugss = 1 - ((alpha*ANC)/(alpha_max*ANC))^gamma;

figure(1);
plot(out.x5);
title("Neutrophils under constant concentration 250mg",'interpreter','latex');
ylabel("Number of neutrophils scaled by $10^9$",'interpreter','latex');
xlabel("Time [days]",'interpreter','latex');
grid on;

figure(2);
plot(out.x5);
title("Neutrophils under constant concentration 176.3mg",'interpreter','latex');
ylabel("Number of neutrophils scaled by $10^9$",'interpreter','latex');
xlabel("Time [days]",'interpreter','latex');
grid on;


%% Part 2: Optimal Dosing

circ0 = 5.04;
ANC = 1.5;
gamma = 0.16;
kcirc = 4/120.4;
ktr = 4/120.4;
kprol = 4/120.4;
alpha = 1;
alpha_max = circ0/ANC;

Edrugi = 0.250;
Edrugf = 1 - ((alpha*ANC)/(alpha_max*ANC))^gamma;
Edrugss = Edrugf;

x5ss = circ0; %circ0*(kprol*(1-Edrugss)/ktr)^(1/gamma);
x4ss = circ0;
x3ss = circ0;
x2ss = circ0;
x1ss = ANC;

Af = [0 0 0 0 -gamma*kprol/ktr;
      1 -1 0 0 0;
      0 1 -1 0 0;
      0 0 1 -1 0;
      0 0 0 1 -kcirc/ktr];
A = ktr*Af;
b1scaled = -kprol*x1ss*((circ0/x5ss)^gamma);
BT = [b1scaled 0 0 0 0];
B = BT';

dxt0 = circ0.*[1;1;1;1;1] - x5ss.*[1;1;1;1;kcirc/ktr];

R2 = 25;
R1 = [100 0 0 0 0;
      0 1 0 0 0;
      0 0 1 0 0;
      0 0 0 1 0;
      0 0 0 0 10];

Fopt = lqr(A,B,R1,R2);

figure(2);
plot(out.time,out.x5);
title("Neutrophils under optimal concentration",'interpreter','latex');
ylabel("Number of neutrophils scaled by $10^9$",'interpreter','latex');
xlabel("Time [days]",'interpreter','latex');
grid on;

fx5 = @(x) (1 - (x/circ0)^gamma);
Edrugst = arrayfun(fx5, out.x5);

figure(3);
plot(out.time,2.5 - 10*Edrugst);
title("Optimal Continuous Drug Schedule",'interpreter','latex');
ylabel("Optimal amount of drug in [mg] scaled by 100",'interpreter','latex');
xlabel("Time [days]",'interpreter','latex');
grid on;
