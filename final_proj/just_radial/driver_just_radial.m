clear; clc;

r_f = 0.00466 ; %[m]
d_g = 0.00003; %[m]
d_c = 0.000673; %[m]
r_c = r_f + d_g + d_c; %[m]
r_g = r_f + d_g; %[m]
k_f = 3 ; %[W/m-K]
k_c = 17; %[W/m-K]
k_g = .25; %[W/m-K]
T_cool = 570; % [K]
h_cool = 25000 ; % [W/m^2*K]
LHR = 20000 ; %[W/m]
Q_max = LHR/(pi*r_f^2) ; % LHR = Q*pi*r_f^2 [W/m^3]

% radial grid
%Number of Grid Points in Each Region
Nf = 15;
Ng = 12;
Nc = 10;

% Radial Heat Generation
% Nuclear Data
SigF = 16.9; %[m^-1]
SigM = 2.2; %[m^-1]
D_f = .0062; %[m]
D_m = 0.00143; %[m]
L_f = sqrt(D_f / SigF); %[m]
L_m = sqrt(D_m / SigM); %[m]
r_m = 0.0063; %[m]

% scratch_work.tex method
S_0 = 1.01916e17;
% S_0 = 3.20179e17; %1.01916e17; %[m^-3 s^-1]
Q_scl = 1.304e-7; %[J/m] 

r_dom = linspace(0, r_f, Nf);
Q_vec = rad_heat_gen(r_dom,Q_scl,r_f,r_m,D_m,D_f,L_f,L_m,S_0);

TCO = T_cool + LHR/(2*pi*h_cool*r_c);

% improved
[T_improved, ~] = radial_solver(Q_vec,LHR,TCO,r_f,r_g,r_c,k_f,k_g,k_c,Nf,Ng,Nc);

% analytical
fuel_grid = linspace(0,r_f,Nf) ;
gap_grid = linspace(r_f,r_g,Ng) ;
clad_grid = linspace(r_g,r_c,Nc) ;

[T_analytical, grid] = fuel_temp_analytical(LHR,r_f,d_g,d_c,k_f,k_g,k_c, ...
    h_cool,T_cool,fuel_grid, gap_grid, clad_grid);
plot(grid,T_analytical,'b',grid,T_improved,'ro')