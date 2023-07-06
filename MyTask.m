%%
clc
clear
close all

c = 0.02;
Sat = 1000;
Ts = 0.02;
%% Define the parameters

mc = 1;
mp = 0.5;
g = 9.82;
l = 0.1;
d = 0.01;

%% Define matrice

A = [0,             0,          1,          0;
     0,             0,          0,          1;
     0,     (g*mp)/mc,     -d/mc, -d/(l*mc);
     0, g*(mc+mp)/(l*mc),-d/(l*mc),   -d*(mc+mp)/(l^2*mc*mp)];
B = [0;
     0;
     1/mc;
     1/(l*mc)];
C = [0,  1,  0,  0;];
D = 0;


  %% Build system
  
sys = ss(A, B, C, D);


%% Discrete Time
sys_d = c2d(sys, Ts);
Ad = sys_d.a;
Bd = sys_d.b;
Cd = sys_d.c;
Dd = sys_d.d;

%% controller
%LQR
Q = 100 * eye(4);
R = 1;

K = lqrd(Ad, Bd, Q, R, Ts);

% P-LQR
K1 = 1e4;
K2 = K1 + 1;


%% full order observer
% des_pole_cd = [0.3; 0.3; 0.3; 0.3];
% L = acker(Ad', Cd',des_pole_cd )';

%% Reduced order observer
P = [Cd;
    0,  0,  1,  0;
    1,  0,  0,  0;
    0,  0,  0,  1];
Ad1 = P*Ad/P;
Bd1 = P*Bd;
Cd1 = C/P;
Ad12 = Ad1([1], [1 2 3]);
Ad22 = Ad1([2 3 4], [2 3 4]);

des_pole_od = (eig(Ad1))/10;
L = acker(Ad22', Ad12', des_pole_od([1, 3, 4]))';






