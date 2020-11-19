function [rad_stress, hoop_stress] = calculate_cladding_stress(r, temp, R_0, R, Nc)
    %%% This function calculates the radial and hoop stress in the fuel
    %%% pellet along a radial profile at constant z.
    
    % Prefactor Values
 
    dr_clad = (R_0-R) / (Nc - 1);
    r_temp = r.*temp;
    int_rT = trapz_integral(r_temp, dr_clad);
    
    rad = zeros(Nc, 1);
    hoop = zeros(Nc, 1);
    
    rad(1) = 0;
    hoop(1) = ( 2 / (R_0^2 - R^2))*int_rT - temp(1);
    
    for i = 2:length(r)
        c_rad = ( 1 - ( R / r(i) )^2 )/( R_0^2 - R^2);
        c_hoop = ( 1 + ( R / r(i) )^2 )/(R_0^2 - R^2);
        int_rTp = (1. / (r(i).^2))* trapz_integral(r_temp(1:i), dr_clad);
        
        rad(i) = c_rad*int_rT - int_rTp;
        hoop(i) = c_hoop*int_rT + int_rTp - temp(i);
    end   
    
    % Coefficient of thermal expansion
    k1 = 6.48E-6 * temp;
    k2 = 1.95E-3*ones(size(temp));
    alf = k1-k2;

    E = 1.088E11 - 5.47E7*temp;
    nu = 0.3;

    units = alf .* (E / (1 - nu));
    rad_stress = units.*rad;
    hoop_stress = units.*hoop;

end