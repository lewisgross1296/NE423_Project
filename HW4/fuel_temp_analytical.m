function T_f = fuel_temp_analytical(r_f,k_f,LHR,T_S,fuel_grid)
    % given the parameters, compute the analytical solution for uniform heat
    % generation temperature distribution in the fuel
    T_f = T_S + LHR/(4*pi*k_f)*(1-(fuel_grid.^2)/r_f^2) ;
end