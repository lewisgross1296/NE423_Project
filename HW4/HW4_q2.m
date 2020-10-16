%% Lewis Gross NE 423 HW 4
clear; clc;

% Question 2

% gauge length in mm
l_0 = 50.8;
% diameter in mm
d = 12.8;
% area in mm squared
A_0 = pi*d^2/4;

loads = [7330 ,15100 ,23100 , 30400 , 34400 , 38400 ,41300 ,44800 , 46200 , ...
          47300 ,47500 , 46100 , 44800 , 42600 , 36400 ];
      
engr_stress = loads / A_0; % N/mm^2, nicely 1 N/mm^2 is 1 MPa
 
% lengths at each point applied load
lengths = [ 50.851, 50.902, 50.952, 51.003 , 51.054 , 51.311 , 51.816 , ...
            52.823 , 53.848 , 54.864 , 55.880 , 56.896 , 57.658 , ...
            58.420 , 59.182 ];
% change in length with applied load
u = lengths - l_0*ones(1,length(lengths));
percent_u = 100* u / l_0;
figure(1);plot(percent_u,engr_stress,'bo',percent_u,engr_stress,'r')
xlabel('Strain \epsilon (%)')
ylabel('Engineering Stress \sigma (N/mm^2)')
title('Engineering Stress-Strain Curve')

% extract elastic region, from inspection, the last fourth point deviattes
% from the line, so we take the first three points to be the elastic region
% and do a linear regression on them and (0,0) 
percent_u_elastic = [0,percent_u(1,1:3)];
engr_stress_elastic = [0, engr_stress(1,1:3)];
figure(2);h=plot(percent_u_elastic,engr_stress_elastic,'c',percent_u_elastic,engr_stress_elastic,'ko');
set(h(1),'linewidth',4) ;
xlabel('Strain \epsilon (%)')
ylabel('Engineering Stress \sigma (N/mm^2)')
title('Elastic Region Linear Plot')

% regression of elastic region to find the modulus of elasticity
% coeffs(1) = slope = modulus of elasticity
coeffs = polyfit(percent_u_elastic,engr_stress_elastic,1);
disp(['The modulus of elasticity is the slope of the elastic region points and is ', ...
       num2str(coeffs(1))])
   
% evaluate that fit a strain offset of 0.2 percent
ys_at_offset = polyval(coeffs, 0.2);
disp(['The yield strength at 0.2 percent offset is ',num2str(ys_at_offset)])

% find the UTS
sig_UTS = max (engr_stress);
disp(['The ultimate tensil strength is ',num2str(sig_UTS)])
