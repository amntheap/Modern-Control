%%
clc
clear
close all

%% Define the parameters

mp = 20;
mw = 6;
l = 0.2;
jpy = 1.33;
r = 0.2;
d = 0.5;
jp = 0.27;
jw = 0.12;
g = 9.82;

%% Define matrice elements-linearized

A23 = -mp^2*l^2*g/(mp*jp+2*(jp+mp*l^2)*(mw+jw/r^2));
A43 = (mp^2*g*l+2*mp*g*l*(mw+jw/r^2))/(mp*jp+2*(jp+mp*l^2)*(mw+jw/r^2));
B21 = (jp+mp*l^2)/(r+mp*l)/(mp*jp+2*(jp+mp*l^2)*(mw+jw/r^2));
B22 = B21;
B41 = (-(r+l)*mp/r-2*(mw+jw/r^2));
B42 = B41;
B61 = d/(2*r)/(jp+(d^2/2*r)*(mw*r+jw/r));
B62 = B61;

%% define matrices

A =[0,  1,  0,  0,  0,  0;
    0,  0,  A23,    0,  0,  0;
    0,  0,  0,  1,  0,  0;
    0,  0,  A43,    0,  0,  0;
    0,  0,  0,  0,  0,  1;
    0,  0,  0,  0,  0,  0];

B = [0,     0;
     B21, B22;
     0,     0;
     B41, B42;
     0,     0;
     B61, B62];

C = [0, 1,  0,  0,  0,  0;
     0, 0,  1,  0,  0,  0;
     0, 0,  0,  0,  0,  1];
 
 D = [0 0;
      0 0;
      0 0];

A1 = [0,     1,      0,       0;
     0,     0,  -6.388,      0;
     0,     0,      0,       1;
     0,     0,   60.6861,    0];
 
 B1 = [0;  0.3803;   0;    -2.3629];
 
 C1 = [0 1 0 0;
       0 0 1 0];
 D1  =[0; 0];
 
 A2 = [0,    1;
      0,    0];
  
  B2 = [0;   0.5085];
  
  C2 = [0 1];
  D2 = 0;
  
  %% Build system
  
  sys = ss(A, B, C, D)
  
  %% Controllability and Observability
  
  Ct = rank(ctrb(sys));
  Ob = rank(obsv(sys));
  
  %% Observable Modes
  [Abaro, Bbaro, Cbaro, To, ko] = obsvf(A, B, C);
  
  %% Controllable Modes
  [Abarc, Bbarc, Cbarc, Tc, kc] = ctrbf(A, B, C);
  
  %% minimal form for the system
  sysr = minreal(sys);
  
  %% State Feedback
  K1 = place(A1, B1, [-2 -3 -4 -5]);
  K2 = place(A2, B2, [-1 -3]);
  %% Far Poles
  %K1 = place(A1, B1, [-20 -30 -40 -50]);
  %K2 = place(A2, B2, [-10 -30]);
  
  %% Robust Tracking
% A_new = [A1 zeros(4,2); C1 zeros(2,2)];
% B_new = [B1;zeros(2,1)];
% C_new = [C1 zeros(2,1)];
% Ki = place(A_new, B_new,[-2 -3 -4 -5 -6 -7]);


%% full order observer
%L1 = place(A1', C1', [-20 -30 -40 -50])';
%L2 = place(A1', C1', [-5 -15])';

%% Reduced order observer

% P = [C1;
%     0,  0,  0,  1;
%     1,  0,  0,  0];
% A_new = P*A1/P;
% B_new = P*B1;
% C_new = C1/P;
% A22_new = A_new([3 4], [3 4]);
% A12_new = A_new([1 2], [3 4]);
% L = place(A22_new', A12_new', [-10 -20])';
%% LQR state feedback
Q1 = eye(4);
R1 = 1;
K1 = lqr(A1, B1, Q1, R1);

Q2 = eye(2);
R2 = 1;
K2 = lqr(A2, B2, Q2, R2);

%% nonlinear model LQR

Q = eye(6);
R = 0.1;
Kn = lqr(A, B, Q, R);
K1 = place(A, B, [-2 -3 -4 -5 -4.5 -3.2]);