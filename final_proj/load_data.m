load('previous_model.mat');

% Scale radial direction to milimeters and set plotting font size
r_scl = r * 1000;
fnt=16;

%%% Identify Z index of of maximumal temperature 
size_temp = size(temp_2D_mesh);
maximum = max(temp_2D_mesh(1, 1:size_temp(2)));
max_idx = find(temp_2D_mesh(1, 1:size_temp(2)) == maximum);

%%% Make temperature color map
figure(2);
[R , Z] = meshgrid(r_scl,z);
surf(R,Z,temp_2D_mesh');
view(2), shading interp
%%% Label Figure
fnt=16;
xlabel('R [mm]','FontSize',fnt);
ylabel('Z [m]','FontSize',fnt);
title(['RZ Temperature Plot. Number of Nodes, fuel = ', num2str(Nf), ', gap = ', num2str(Ng), ', cladding = ', num2str(Nc), ', Z points = ', num2str(M)], 'FontSize',fnt+2);
a = colorbar('EastOutside');
a.Label.String = 'Temperature [K]';
a.Label.FontSize = fnt;



