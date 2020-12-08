load('previous_model.mat');
temp_previous = temp_2D_mesh;
load('improved_model.mat');
temp_improved = temp_2D_mesh;

delta = temp_improved - temp_previous;

r_scl = r * 1000;
[R , Z] = meshgrid(r_scl,z);
surf(R,Z,delta');
view(2), shading interp
%%% Label Figure
fnt=16;
xlabel('R [mm]','FontSize',fnt);
ylabel('Z [m]','FontSize',fnt);
title(['Difference in Models. Number of Nodes, fuel = ', num2str(Nf), ', gap = ', num2str(Ng), ', cladding = ', num2str(Nc), ', Z points = ', num2str(M)], 'FontSize',fnt+2);
a = colorbar('EastOutside');
a.Label.String = 'Temperature [K]';
a.Label.FontSize = fnt;