%% Lewis Gross NE 423 HW 4

% Question 3
R = 0.535; % cm
R0 = 0.605; % cm
rat = R0/R ;
T0 = 293 ;% K   
Tstop = 1628 ; % K guess past break
Tvec = [T0:5:Tstop] ; % 5 K per second
times  = [0:1:(Tstop - T0)/5] ;% change in temp / heating rate = final time
if(length(Tvec)~=length(times))
    disp('issue with time/temp mapping' )
end
sig_UTS = 310 - 0.17*Tvec ;

% starting at 1 MPa
Pvec_1 = Tvec.*(1/T0); % MPa
sig_equiv1 = Pvec_1 * sqrt(rat^2+(rat+1)^2+1) / ( 2*sqrt(2)*(rat-1) );
figure(1);plot(Tvec,sig_equiv1-sig_UTS,'bo',Tvec,zeros(size(Tvec)),'k')
xlabel('Temperature [K]')
ylabel('\sigma_{equiv} - \sigma_{UTS} [MPa]')
title('Comparing Equivalent Stress and UTS to Find Failure Temperature,P0=1 MPa')


% starting at 4 MPa
Pvec_4 = Tvec.*(4/T0); % MPa
sig_equiv4 = Pvec_4 * sqrt(rat^2+(rat+1)^2+1) / ( 2*sqrt(2)*(rat-1) );
figure(2);plot(Tvec,sig_equiv4-sig_UTS,'ro',Tvec,zeros(size(Tvec)),'k')
xlabel('Temperature [K]')
ylabel('\sigma_{equiv} - \sigma_{UTS} [MPa]')
title('Comparing Equivalent Stress and UTS to Find Failure Temperature,P0=4 MPa')

% starting at 7 MPa
Pvec_7 = Tvec.*(7/T0); % MPa
sig_equiv7 = Pvec_7 * sqrt(rat^2+(rat+1)^2+1) / ( 2*sqrt(2)*(rat-1) );
figure(3);plot(Tvec,sig_equiv7-sig_UTS,'go',Tvec,zeros(size(Tvec)),'k')
xlabel('Temperature [K]')
ylabel('\sigma_{equiv} - \sigma_{UTS} [MPa]')
title('Comparing Equivalent Stress and UTS to Find Failure Temperature,P0=7 MPa')
