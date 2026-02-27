clc; clear; close all;

%% Input Signal

f0 = 1000;
fs = 96000;
Ts = 1/fs;
StopTime = 1;
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
Z_D1 = 1;
Z_M1 = 1;
Z_M2 = 1;

Z_dist = diag([Z_C5, R11, Rin, Z_D1, R9, Rin, Z_M2, R10, R7, Z_C6, R8, Z_M1, Z_C4]);

Bv_dist = [-1 0  0  0  0  0  0  1  0  0  0  0  0;
           0  0  1  0  0 -1  1  0  1  0  0  0  0;
           0  1  0  0 -1 -1  1  0  0  1  0  0  0;
           0  0  0  0 -1  0  1  0  0  0  1  0  0;
           0  0  0  1  0  0  1  0  0  0  0  1  0;
           0  0  0  0  0  0 -1  0  0  0  0  0  1];

Bi_dist = [-1  0  0  0  0  0  0  1  0  0  0  0  0;
            0  0  1  0  0 -1  1  0  1  0  0  0  0;
            1  1  0  0  0  0  0  0  0  1  0  0  0;
            0  0  0  0 -1  0  1  0  0  0  1  0  0;
            0  0  0  1  0  0  1  0  0  0  0  1  0;
            0  0  0  0  0  0 -1  0  0  0  0  0  1];

S_dist = eye(13) - 2*Z_dist*Bi_dist'*((Bv_dist*Z_dist*Bi_dist')\Bv_dist);

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

% Both sorted in increasing order
load("diode_char.mat");
load("mos_char.mat");

slope_d = (Vvect_d(2:end) - Vvect_d(1:end-1) + 1e-12)./(Ivect_d(2:end) - Ivect_d(1:end-1) + 1e-12)  + 1e-12;
slope_mos = (Vvect(2:end) - Vvect(1:end-1) + 1e-12)./(Ivect(2:end) - Ivect(1:end-1) + 1e-12) + 1e-12;

% Kirchhoff to Wave conversion
aVect_m1 = (Vvect + Z_dist(12,12) * Ivect);
bVect_m1 = (Vvect - Z_dist(12,12) * Ivect);

aVect_m2 = (Vvect + Z_dist(7,7) * Ivect);
bVect_m2 = (Vvect - Z_dist(7,7) * Ivect);

aVect_d = (Vvect_d + Z_dist(4,4) * Ivect_d);
bVect_d = (Vvect_d - Z_dist(4,4) * Ivect_d);

[mu0_m1, mu1_m1, etaVect_m1] = CPWL_param(aVect_m1, bVect_m1, 0); 
[mu0_m2, mu1_m2, etaVect_m2] = CPWL_param(aVect_m2, bVect_m2, 0);
[mu0_d, mu1_d, etaVect_d] = CPWL_param(aVect_d, bVect_d, 0);

%% Algorithm

% Initialization
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
    b_dist(1) = (a_dist(1) + b_dist(1))/2;
    b_dist(3) = V_1(n);
    b_dist(6) = V_ref;
    b_dist(10) = (a_dist(10) + b_dist(10))/2;
    b_dist(13) = (a_dist(13) + b_dist(13))/2;
    
    % Nonlinearties
    flag = 1;
    S_temp = S_dist;

    while flag

        % (Nonlinear) Local Scattering Stage
        
        [b_dist(4), idx3]= CPWL_function(a_dist(4), aVect_d, mu0_d, mu1_d, etaVect_d);
        [b_dist(7), idx2]= CPWL_function(a_dist(7), aVect_m2, mu0_m2, mu1_m2, etaVect_m2);
        [b_dist(12), idx1]= CPWL_function(a_dist(12), aVect_m1, mu0_m1, mu1_m1, etaVect_m1);
        
        % Global Scattering Stage 

        a_dist = S_temp * b_dist; 

        % Convergence Check

        v = (a_dist + b_dist)/2;

        if max(abs(v - v_old_iter)) < eSIM

            % Update of Port Resistances
            Z_dist(4, 4) = slope_d(idx3);
            Z_dist(7, 7) = slope_mos(idx2);
            Z_dist(12, 12) = slope_mos(idx1);

            aVect_m1 = (Vvect + Z_dist(12,12) * Ivect);
            bVect_m1 = (Vvect - Z_dist(12,12) * Ivect);

            aVect_m2 = (Vvect + Z_dist(7,7) * Ivect);
            bVect_m2 = (Vvect - Z_dist(7,7) * Ivect);

            aVect_d = (Vvect_d + Z_dist(4,4) * Ivect_d);
            bVect_d = (Vvect_d - Z_dist(4,4) * Ivect_d);

            [mu0_m1, mu1_m1, etaVect_m1] = CPWL_param(aVect_m1, bVect_m1, 0); 
            [mu0_m2, mu1_m2, etaVect_m2] = CPWL_param(aVect_m2, bVect_m2, 0);
            [mu0_d, mu1_d, etaVect_d] = CPWL_param(aVect_d, bVect_d, 0);

            flag = 0;     
        end  

        v_old_iter = v;
    end

    a_dist = S_dist * b_dist;

    S_dist = eye(13) - 2*Z_dist*Bi_dist'*((Bv_dist*Z_dist*Bi_dist')\Bv_dist);

    V_2(n) = (a_dist(8) + b_dist(8))/2 + (a_dist(10) + b_dist(10))/2 + (a_dist(2) + b_dist(2))/2;

    % Output Stage
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

