dx = 0.01;
x = 0:dx:.5*pi;
F = cos(x);

plot(x, F)
trapz_integral(F, dx)