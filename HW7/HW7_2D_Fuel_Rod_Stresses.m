%% NE 423 HW5 ~Computational Part~
% Aysia Demby
% Michael Gerard
% Lewis Gross

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
Nf = 65;
Ng = 35;
Nc = 35;

% axial grid
H = 4; % m
% Number of Grid Points in the z direction, use odd number
M = 65;
dz = H./(M-1);
mdot = 0.25 ;% kg/s
CPW = 4200 ;% J/kg-K

N = Nf+Ng+Nc-2 ; 
temp_2D_mesh = zeros(N,M);
z = [-H/2:dz:H/2]';
% compute the radial temperature distribution for each z
for j = 1:length(z)
    Q = vol_heat_gen(Q_max,z(j),H);
    TCO = cladding_outer(z(j),pi*r_f^2,Q,r_c,mdot,CPW,h_cool,T_cool,H);
    [T , r] = radial_solver(Q,TCO,r_f,r_g,r_c,k_f,k_g,k_c,Nf,Ng,Nc);
    temp_2D_mesh(:,j) = T
end

save('fine_mesh_HW7')