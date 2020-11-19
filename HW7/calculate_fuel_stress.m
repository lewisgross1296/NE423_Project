function [rad_stress, hoop_stress] = calculate_fuel_stress(r, temp, r_f, Nf)
    %%% This function calculates the radial and hoop stress in the fuel
    %%% pellet along a radial profile at constant z.
    
    % Constant Values
    rf_sqrd = 1. / (r_f * r_f);
    
    dr_fuel = r_f / (Nf - 1);
    r_temp = r.*temp;
    int_rT = rf_sqrd * trapz_integral(r_temp, dr_fuel);
    
    rad = zeros(Nf, 1);
    hoop = zeros(Nf, 1);
    
    rad(1) = int_rT;
    hoop(1) = int_rT - temp(1);
    
    for i = 2:length(r)
 
        int_rTp = (1. / (r(i).^2)) * trapz_integral(r_temp(1:i), dr_fuel);
        
        rad(i) = int_rT - int_rTp;
        hoop(i) = int_rT + int_rTp - temp(i);
    end
    
    % Coefficient of thermal expansion
    k1 = 1e-5 * temp;
    k2 = 3e-3*ones(size(temp));
    k3 = zeros(size(temp));
    for i = 1:Nf
        k3(i) = 0.04 * exp( -6.9e-20 / (1.3807e-23 * temp(i)) );
    end
    alf = k1 - k2 + k3;

    E = 2.3e11;
    nu = 0.316;

    units = alf * (E / (1 - nu));
    rad_stress = units.*rad;
    hoop_stress = units.*hoop;

end