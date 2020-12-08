function [T_vec, grid] = fuel_temp_analytical(LHR,r_f,d_g,d_c,k_f,k_g, ...
                         k_c,h_cool,T_cool,fuel_grid, gap_grid, clad_grid)
    % given the parameters, compute the analytical solution for uniform heat
    % generation temperature distribution in the fuel
    % radii
    r_g = r_f + d_g; % gap outer radius
    r_c = r_g + d_c; % cladding outer radius

    T_CO = T_cool + LHR/(2*pi*h_cool*r_c);
    T_CI = T_CO + LHR*log(1+(d_c/r_f))/(2*pi*k_c);
    T_S = T_CI + LHR*d_g/(2*pi*k_g*r_f);

    T_f = T_S + LHR/(4*pi*k_f)*(1-(fuel_grid.^2)/r_f^2)  ;
    T_g = T_S - LHR/(2*pi*k_g)*log(gap_grid/r_f)  ;
    T_c = T_CI - LHR/(2*pi*k_c)*log(clad_grid/(r_f+d_g)) ;
    T_vec = [T_f(1:end-1) T_g T_c(2:end)]' ;
    grid  =  [fuel_grid(1:end-1) gap_grid clad_grid(2:end)]';
end