function TCO = cladding_outer(z,Af,Qmax,rc,mdot,C,h,Tcool_in,H)
    % z,Af,Qmax,rc,mdot,C,h,Tcool_in TODO
    TCO = Tcool_in + (Af*Qmax)./(pi*mdot*C).*(sin(pi.*z./H)+1) + (Af*Qmax)./(h*2*pi.*rc).*cos(pi.*z./H);
end