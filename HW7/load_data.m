load('fine_mesh_HW7');

% Scale radial direction to milimeters and set plotting font size
r_scl = r * 1000;
fnt=16;

%%% Identify Z index of of maximumal temperature 
size_temp = size(temp_2D_mesh);
maximum = max(temp_2D_mesh(1, 1:size_temp(2)));
max_idx = find(temp_2D_mesh(1, 1:size_temp(2)) == maximum);

%%% Determine shifted indices along Z axis
raise_idx = max_idx + round(0.8 * (size_temp(2) - max_idx));
lower_idx = max_idx - round(0.4 * (size_temp(2) - max_idx));

%%% Calculate Fuel Stresses
[fuel_rad_stress_1, fuel_hoop_stress_1] = calculate_fuel_stress(r(1:Nf), temp_2D_mesh(1:Nf, max_idx), r_f, Nf);
[fuel_rad_stress_2, fuel_hoop_stress_2] = calculate_fuel_stress(r(1:Nf), temp_2D_mesh(1:Nf, raise_idx), r_f, Nf);
[fuel_rad_stress_3, fuel_hoop_stress_3] = calculate_fuel_stress(r(1:Nf), temp_2D_mesh(1:Nf, lower_idx), r_f, Nf);

% Calculate Cladding Stresses
[clad_rad_stress_1, clad_hoop_stress_1] = calculate_cladding_stress(r(Nf + Nc - 1:N), temp_2D_mesh((Nf + Nc - 1:N), max_idx), r_c, r_g, Nc);
[clad_rad_stress_2, clad_hoop_stress_2] = calculate_cladding_stress(r(Nf + Nc - 1:N), temp_2D_mesh((Nf + Nc - 1:N), raise_idx), r_c, r_g, Nc);
[clad_rad_stress_3, clad_hoop_stress_3] = calculate_cladding_stress(r(Nf + Nc - 1:N), temp_2D_mesh((Nf + Nc - 1:N), lower_idx), r_c, r_g, Nc);

%%% Plot Fuel Radial Stress
figure(1);
plot(r_scl(1:Nf), fuel_rad_stress_1*1e-9,'color','k','linewidth',3);
hold on
plot(r_scl(1:Nf), fuel_rad_stress_2*1e-9,'color','r','linewidth',3);
plot(r_scl(1:Nf), fuel_rad_stress_3*1e-9,'color','g','linewidth',3);
hold off

grid
xlabel('R [mm]','FontSize',fnt);
ylabel('$\sigma_{\theta}^{th}$ [GPa]','Interpreter','latex','FontSize',fnt);
title('Radial Stress in Fuel','FontSize',fnt+2);
legend('max temp stress','lower half','upper half')
figure;

%Plot Cladding Radial Stress
figure(2);
plot(r_scl(Nf + Nc - 1:N), clad_rad_stress_1*1e-9,'color','k','linewidth',3);
hold on
plot(r_scl(Nf + Nc - 1:N), clad_rad_stress_2*1e-9,'color','r','linewidth',3);
plot(r_scl(Nf + Nc - 1:N), clad_rad_stress_3*1e-9,'color','g','linewidth',3);
hold off

grid
xlabel('R [mm]','FontSize',fnt);
ylabel('$\sigma_{r}^{th}$ [GPa]','Interpreter','latex','FontSize',fnt);
title('Radial Stress in Cladding','FontSize',fnt+2);
legend('max temp stress','lower half','upper half')
figure;

%%% Plot Fuel Hoop Stress
figure(3);
plot(r_scl(1:Nf), fuel_hoop_stress_1*1e-9,'color','k','linewidth',3);
hold on
plot(r_scl(1:Nf), fuel_hoop_stress_2*1e-9,'color','r','linewidth',3);
plot(r_scl(1:Nf), fuel_hoop_stress_3*1e-9,'color','g','linewidth',3);
hold off

grid
xlabel('R [mm]','FontSize',fnt);
ylabel('$\sigma_{\theta}^{th}$ [GPa]','Interpreter','latex','FontSize',fnt);
title('Hoop Stress in Fuel','FontSize',fnt+2);
legend('max temp stress','lower half','upper half')
figure;

%%Plot Cladding Hoop Stress
figure(4);
plot(r_scl(Nf + Nc - 1:N), clad_hoop_stress_1*1e-9,'color','k','linewidth',3);
hold on
plot(r_scl(Nf + Nc - 1:N), clad_hoop_stress_2*1e-9,'color','r','linewidth',3);
plot(r_scl(Nf + Nc - 1:N), clad_hoop_stress_3*1e-9,'color','g','linewidth',3);
hold off

grid
xlabel('R [mm]','FontSize',fnt);
ylabel('$\sigma_{\theta}^{th}$ [GPa]','Interpreter','latex','FontSize',fnt);
title('Hoop Stress in Cladding','FontSize',fnt+2);
legend('max temp stress','lower half','upper half')
figure;

%%% Make temperature color map
figure(5);
[R , Z] = meshgrid(r_scl,z);
surf(R,Z,temp_2D_mesh');
view(2), shading interp;

hold on

%%% Make Z slice arrays
z_val_1 = z(max_idx) * ones(N, 1);
z_val_2 = z(raise_idx) * ones(N, 1);
z_val_3 = z(lower_idx) * ones(N, 1);

%%% Overlay Z slice arrays onto temerature color map
z_order = (max(temp_2D_mesh(1,1:size_temp(2)))+1) * ones(N, 1);

plot3(r_scl, z_val_1, z_order, 'color','k','linewidth',3)
plot3(r_scl, z_val_2, z_order, 'color','r','linewidth',3)
plot3(r_scl, z_val_3, z_order, 'color','g','linewidth',3)

%%% Label Figure
fnt=16;

xlabel('R [mm]','FontSize',fnt);
ylabel('Z [m]','FontSize',fnt);

title(['RZ Temperature Plot. Number of Nodes, fuel = ', num2str(Nf), ', gap = ', num2str(Ng), ', cladding = ', num2str(Nc), ', Z points = ', num2str(M)], 'FontSize',fnt+2);
a = colorbar('EastOutside');
a.Label.String = 'Temperature [K]';
a.Label.FontSize = fnt;

saveas(figure(1),'Fuel_Radial_Stress.png');
saveas(figure(2),'Clad_Radial_Stress.png');
saveas(figure(3),'Fuel_Hoop_Stress.png');
saveas(figure(4),'Clad_Hoop_Stress.png');
saveas(figure(5),'Temp_Plot.png');