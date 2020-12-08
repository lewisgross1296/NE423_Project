function Q = rad_heat_gen(r, Qscl, Rf, Rm, Dm, Df, Lf, Lm, S0)
    Rf_Lf = Rf / Lf;
    Rf_Lm = Rf / Lm;
    Rm_Lm = Rm / Lm;
    
    numer = ( besseli(0,Rf_Lm) / besseli(0,Rf_Lf) ) +...
        ( ( besseli(1,Rm_Lm) * besselk(0,Rf_Lm) ) / ( besselk(1,Rm_Lm) * besseli(0,Rf_Lf) ) );
    denom =  ((Dm * Lf) / (Df * Lm)) * ( ( besseli(1, Rf_Lm) / besseli(1, Rf_Lf) ) -...
        ( ( besseli(1, Rm_Lm) * besselk(1, Rf_Lm) ) / ( besselk(1, Rm_Lm) * besseli(1, Rf_Lf) ) ) ) -...
        ( ( besseli(0, Rf_Lm) / besseli(0, Rf_Lf) ) + ( ( besseli(1, Rm_Lm) * besselk(0, Rf_Lm) ) / ( besselk(1, Rm_Lm) * besseli(0, Rf_Lf) ) ) );
    zeta = ((Lm * Lm) / (Dm * besseli(0,Rf_Lf))) * ((numer / denom) + 1);
    
    phi = zeta * S0 * besseli(0, r / Lf);
    Q = Qscl * phi;
end