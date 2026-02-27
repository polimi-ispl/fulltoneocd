clc; clear; close all;

%% Input Signal

f0 = 1000;
fs = 96000;
Ts = 1/fs;
StopTime = 5;
t = 0:Ts:StopTime-Ts;
Vin = sin(2*pi*f0*t);
N = length(Vin);

V_ref = 4.5;

%% Input Stage

Rin = 1e-6;
R1 = 10e3;
R2 = 1e6;
R3 = 470e3;
R4 = 18e3;
R5 = 500e3;
R6 = 2.2e3;
Z_C1 = Ts/(22e-9); % BE
Z_C2 = Ts/(220e-12); % BE
Z_C3 = Ts/(100e-9); % BE

Z_in = diag([R5, Z_C2, R6, Rin, R1, Z_C1, R2, R3, R4, Z_C3]);

% Fundamental Loop Matrix of the V-net
Bv = [0  0  0  1  0  0  1  0  0  0;
      0  0  0  1  1  1  0  1  0  0;
      1 -1  0  0  0  0  0  0  1  0;
      0  0  1 -1 -1 -1  0  0  0  1];

% Fundamental Loop Matrix of the I-net
Bi = [ 0  0  0  1  0  0  1  0  0  0;
       0  0  0  1  1  1  0  1  0  0;
       1 -1  0  0  0  0  0  0  1  0;
       0  1  1  0  0  0  0  0  0  1];

S_in = eye(10) - 2*Z_in*Bi'*((Bv*Z_in*Bi')\Bv);

%% Clipping Stage

R7 = 10e3; 
R8 = 220e3;
R9 = 10e3;
R11 = 39e3;
R10 = 150e3;
Z_C4 = Ts/(10e-9); % BE
Z_C5 = Ts/(220e-12); % BE
Z_C6 = Ts/(100e-9); % BE

switch fs
    case 48000
        Z_NL = 1.711309523868095e+03; % from adaptation condition (fs=48k and BE)
    case 96000    
        Z_NL = 9.395424836777855e+02; % from adaptation condition (fs=96k and BE)
end

Z_dist = diag([R10, R11, R7, R9, R8, Rin, Z_C5, Rin, Z_C6, Z_NL, Z_C4]);

Bv_dist = [-1 0  0  0  0  0  1  0  0  0  0;
           0  0  1  1 -1 -1  0  1  0  0  0;
           0  1  0  0 -1 -1  0  0  1  0  0;
           0  0  0 -1  1  0  0  0  0  1  0;
           0  0  0 -1  1  0  0  0  0  0  1];

Bi_dist = [-1  0  0  0  0  0  1  0  0  0  0;
            0  0  1  1 -1 -1  0  1  0  0  0;
            1  1  0  0  0  0  0  0  1  0  0;
            0  0  0 -1  1  0  0  0  0  1  0;
            0  0  0 -1  1  0  0  0  0  0  1];

S_dist = eye(11) - 2*Z_dist*Bi_dist'*((Bv_dist*Z_dist*Bi_dist')\Bv_dist);

%% Output Stage

R14 = 10e3;
R15 = 500e3;
R12 = 22e3;
R13 = 33e3;
Z_C7 = Ts/(10e-6); % BE
Z_C8 = Ts/(47e-9); % BE


Z_out = diag([R14, R12, Rin, Z_C8, R15, Z_C7, R13]);

B_out = [ 1  0  0  1  1  0  0;
          1  1  1  1  0  1  0;
          0 -1  0  0  0  0  1];

S_out = eye(7) - 2*Z_out*B_out'*((B_out*Z_out*B_out')\B_out);

%% Initialization of CPWL Functions

% Sorted in increasing order
load("singleNL_char.mat");

% Kirchhoff to Wave conversion
aVect = (Vvect_s + Z_dist(10, 10) * Ivect_s);
bVect = (Vvect_s - Z_dist(10, 10) * Ivect_s);

[mu0, mu1, etaVect] = CPWL_param(aVect, bVect, 0); 

%% Algorithm

% Initialization of arrays
b_in = zeros(length(Z_in), 1);
a_in = zeros(length(Z_in), 1);
V_1 = zeros(N, 1);

b_dist = zeros(length(Z_dist), 1);
a_dist = zeros(length(Z_dist), 1);
V_2 = zeros(N, 1);

b_out = zeros(length(Z_out), 1);
a_out = zeros(length(Z_out), 1);
V_out = zeros(N, 1);

v_old_iter = 0;
eSIM = 1e-6;

for n = 1:N
    
    % Input stage
    b_in(2) = (a_in(2) + b_in(2))/2;
    b_in(4) = Vin(n);
    b_in(6) = (a_in(6) + b_in(6))/2;
    b_in(8) = V_ref;
    b_in(10) = (a_in(10) + b_in(10))/2;

    a_in = S_in * b_in;

    V_1(n) = (a_in(2) + b_in(2))/2 + (a_in(3) + b_in(3))/2 - (a_in(10) + b_in(10))/2;

    % Distortion stage
    b_dist(6) = V_ref;
    b_dist(7) = (a_dist(7) + b_dist(7))/2;
    b_dist(8) = V_1(n);
    b_dist(9) = (a_dist(9) + b_dist(9))/2;
    b_dist(11) = (a_dist(11) + b_dist(11))/2;

    a_dist(10) = S_dist(10, :) * b_dist;

    % Root Scattering Stage
    [b_dist(10), idx1]= CPWL_function(a_dist(10), aVect, mu0, mu1, etaVect);

    a_dist = S_dist * b_dist;

    V_2(n) = (a_dist(1) + b_dist(1))/2 + (a_dist(9) + b_dist(9))/2 + (a_dist(2) + b_dist(2))/2;

    % Output stage
    b_out(3) = V_2(n);
    b_out(4) = (a_out(4) + b_out(4))/2;
    b_out(6) = (a_out(6) + b_out(6))/2;
    
    a_out = S_out * b_out;
    
    V_out(n) = (a_out(5) + b_out(5))/2;
    
end

%% Plot

figure('Color', 'white')
plot(t, V_out, 'b','LineWidth',2,'DisplayName','WD');
xlim([0.995, 1]);
xlabel('Time [s]','interpreter','latex','FontSize',18);
ylabel('Voltage [V]','interpreter','latex','FontSize',18);
legend('show','interpreter','latex','FontSize',13);
