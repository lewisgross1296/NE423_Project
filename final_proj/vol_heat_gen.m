function Q = vol_heat_gen(Qmax,z,H)
    % accounts for z dependence of flux to scale the Q(r)  
    Q = Qmax*cos(z.*pi./H);
end