function I = trapz_integral(F, dx)
    % This function perfoms a numerical integration on the input F array
    % using the trapezoidal rule.
    twos = 2.*ones(1, length(F)-2);
    weights = [1 twos 1];
    I = 0.5 * dx * sum( dot(F, weights.') );
end