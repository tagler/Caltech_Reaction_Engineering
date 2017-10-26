
% ------------------------------------------------------------------------
% REACTION ENGINEERING FINAL PROJECT
% ------------------------------------------------------------------------
 
% ------------------------------------------------------------------------
% Cocurrent System (Mode = 1)
% ------------------------------------------------------------------------

clear all

% Initial Conditions
X1o = 0.1;
X2o = 0.0;
X3o = 0.07;
X4o = 0.0;
X5o = 0.02;
TRo = 353; % K
TCo = 350; % K
MWo = 37.9; % lb/lbmol
PTo = 114.7; % psia

% Z Range and Initial Condition Vectors
zspan = 0:0.25:20;
Y0 = [X1o; X2o; X3o; X4o; X5o; TRo; TCo; MWo; PTo];
YP0 = [0; 0; 0; 0; 0; 0; 0; 0; 0];

% Preallocate Tmatrix to optimize loop
ss = 0.2; % degrees K
msize = 59;
Tmatrix=ones(msize,2);

for n=1:msize % loop for various TCo temps (350-365)

    % Solve System
    [Y0,YP0] = decic('systemeqs',0,Y0,ones(9,1),YP0,zeros(9,1));
    options = odeset('RelTol',1e-2);
    [z,y] = ode15i('systemeqs',zspan,Y0,YP0,options);
    Tmatrix(n,1) = Y0(7);
    Tmatrix(n,2) = max(y(:,6));
    Y0(7) = Y0(7) + ss;

end

% Plot X1-X5
figure
plot(z,y(:,1),'b.-');
hold
plot(z,y(:,2),'r.-');
plot(z,y(:,3),'g.-');
plot(z,y(:,4),'k.-');
plot(z,y(:,5),'c.-');
title('Cocurrent Composition vs. Position');
ylabel('X');
xlabel('Z (ft)')
legend('DT','VA','O2','CO2','H20');

% Plot TC and TR
figure
plot(z,y(:,6),'r.-');
hold
plot(z,y(:,7),'b.-');
title('Cocurrent Temperature vs. Position');
ylabel('T (K)');
xlabel('Z (ft)')
legend('TR','TC');

% Plot MW
figure
plot(z,y(:,8),'b.-');
title('Cocurrent Mean Molecular Weight vs. Position');
ylabel('MW (lb/lbmol)');
xlabel('Z (ft)')

% Plot PT
figure
plot(z,y(:,9),'b.-');
title('Cocurrent Total Pressure vs. Position');
ylabel('P (psia)');
xlabel('Z (ft)')


% ------------------------------------------------------------------------
% Counter Current System (mode = -1)
% ------------------------------------------------------------------------

% Initial Conditions
X1o = 0.1;
X2o = 0.0;
X3o = 0.07;
X4o = 0.0;
X5o = 0.02;
TRo = 353; % K
TCo = 352.75; % K
MWo = 37.9; % lb/lbmol
PTo = 114.7; % psia

% Z Range and Initial Condition Vectors
zspan = 0:0.25:20;
Y0 = [X1o; X2o; X3o; X4o; X5o; TRo; TCo; MWo; PTo];
YP0 = [0; 0; 0; 0; 0; 0; 0; 0; 0];

% Preallocate Tmatrix to optimize loop
ss = 0.2;
msize = 62 ;
Tmatrix2=ones(msize,2);

for n=1:msize % loop for various TCo temps (350-365)

    % Solve System
    [Y0,YP0] = decic('systemeqs2',0,Y0,ones(9,1),YP0,zeros(9,1));
    options = odeset('RelTol',1e-2);
    [z,y] = ode15i('systemeqs2',zspan,Y0,YP0,options);
    Tmatrix2(n,1) = min(y(:,7));
    Tmatrix2(n,2) = max(y(:,6));
    Y0(7) = Y0(7) + ss;

end

% Plot X1-X5
figure
plot(z,y(:,1),'b.-');
hold
plot(z,y(:,2),'r.-');
plot(z,y(:,3),'g.-');
plot(z,y(:,4),'k.-');
plot(z,y(:,5),'c.-');
title('Countercurrent Composition vs. Position');
ylabel('X');
xlabel('Z (ft)')
legend('DT','VA','O2','CO2','H20');

% Plot TC and TR
figure
plot(z,y(:,6),'r.-');
hold
plot(z,y(:,7),'b.-');
title('Countercurrent Temperature vs. Position');
ylabel('T (K)');
xlabel('Z (ft)')
legend('TR','TC');

% Plot MW
figure
plot(z,y(:,8),'b.-');
title('Countercurrent Mean Molecular Weight vs. Position');
ylabel('MW (lb/lbmol)');
xlabel('Z (ft)')

% Plot PT
figure
plot(z,y(:,9),'b.-');
title('Countercurrent Total Pressure vs. Position');
ylabel('P (psia)');
xlabel('Z (ft)')

% Plot Coolant vs. Fluid Temp
figure
plot( Tmatrix(:,1), Tmatrix(:,2), 'b.-');
hold
plot( Tmatrix2(:,1), Tmatrix2(:,2), 'r.-');
xlabel('Inlet Coolant Temperature, TCo (K)');
ylabel('Maximum Reacting Fluid Temperature, TR (K)');
title('Inlet Coolant Temperature vs. Maximum Reacting Fluid Temperature');
legend('Cocurrent Flow','Countercurrent Flow');
function dy = systemeqs(z,y,yp)

% -------------------------------------------------------------------------
% Cocurrent System (mode = 1) ---------------------------------------------
% -------------------------------------------------------------------------

% Define f(z,y,yp)=0 vector, implicit diff eq system
dy = ones(9,1);

% -------------------------------------------------------------------------
% Rates (lbmol / lbcat h) -------------------------------------------------
% -------------------------------------------------------------------------

% RT and Partial Pressure Expressions
RT = ( (1/y(6))-(1/373) )/1.987;
PDT = y(1)*y(9); %psia
PVA = y(2)*y(9); %psia
PO2 = y(3)*y(9); %psia

% Reaction Constants (k1-k12)
k = [ 1.771*10^-3; 23295.0; 0.5; 1.0; 0.8184; 0.5; 0.2314;
    1.0; 1.25; 2.0; 2.795*10^-4; 33000.0; 0.5; 2.0];

% Reaction 1 Rate Expression
r1top = (k(1)*exp(-k(2)*RT)*(PO2^k(3))*(PDT^k(4)));
r1bom =(1 + k(5)*(PO2^k(6))+ k(7)*(PDT^k(8))+ k(9)*(PVA^k(10)));
ra = r1top / r1bom; % lbmol/(lbcat-h)

% Reaction 2 Rate Expresion
r2top = (k(11)*exp(-k(12)*RT)*(PO2^k(13))*(PVA^k(14)));
r2bom =(1 + k(5)*(PO2^k(6))+ k(7)*(PDT^k(8))+ k(9)*(PVA^k(10)));
rb = r2top / r2bom; % lbmol/(lbcat-h)

R1 = -ra;
R2 = ra - rb;
R3 = -0.5*ra - 2.5*rb;
R4 = 2*rb;
R5 = ra + 2*rb;

% -------------------------------------------------------------------------
% Mass Balances: dy(1)-dy(5)  ---------------------------------------------
% -------------------------------------------------------------------------

% Constants
NT = 2500; % number tubes
rhoB = 100; %lb/ft3
dt = 1/12; % ft
rt = dt/2; % ft
Ac = pi*(rt^2); %ft^2
W = 100000; % lb/h

% dy(1) = dX1, X1-X5 unitless
dy1f1 = ((NT*rhoB*Ac*y(8))/W)*(R1);
dy1f2 = (y(1)/y(8))*yp(8);
dy(1) = dy1f1 + dy1f2 - yp(1);

% dy(2) = dX2, X1-X5 unitless
dy2f1 = ((NT*rhoB*Ac*y(8))/W)*(R2);
dy2f2 = (y(2)/y(8))*yp(8);
dy(2) = dy2f1 + dy2f2 - yp(2);

% dy(3) = dX3, X1-X5 unitless
dy3f1 = ((NT*rhoB*Ac*y(8))/W)*(R3);
dy3f2 = (y(3)/y(8))*yp(8);
dy(3) = dy3f1 + dy3f2 - yp(3);

% dy(4) = dX4, X1-X5 unitless
dy4f1 = ((NT*rhoB*Ac*y(8))/W)*(R4);
dy4f2 = (y(4)/y(8))*yp(8);
dy(4) = dy4f1 + dy4f2 - yp(4);

% dy(5) = dX5, X1-X5 unitless
dy5f1 = ((NT*rhoB*Ac*y(8))/W)*(R5);
dy5f2 = (y(5)/y(8))*yp(8);
dy(5) = dy5f1 + dy5f2 - yp(5);

% -------------------------------------------------------------------------
% Energy Balances: dy(6)-dy(7) --------------------------------------------
% -------------------------------------------------------------------------

% Constants
mode = 1; % type of flow
U = 120.5; % BTU/(h-ft^2-C)
FcCpc = 10^6; % BTU/(h C)
Hr1 = -74.32*1000; % BTU/lbmol
Hr2 = -474.57*1000; % BTU/lbmol
Cp = 0.50; % BTU/(lb C)

% dy(6) = dTR, TR (K)
dy6f1 = (-1*NT)*Ac*(Hr1*ra + Hr2*rb)*rhoB;
dy6f2 = pi*NT*dt*U*( y(6)-y(7) );
dy6f3 = W*Cp;
dy(6) = ((dy6f1 - dy6f2)/dy6f3) - yp(6);

% dy(7) = dTC, TC (K)
dy(7) = ((mode)*pi*NT*dt*U*(y(6)-y(7))) / ((FcCpc)) - yp(7);

% -------------------------------------------------------------------------
% Molecular Weight: dy(8) -------------------------------------------------
% -------------------------------------------------------------------------

% Molecular Weights
M = [46; 44; 32; 44; 18]; % lb/lbmol

% dy(8) = dMW (lb/lbmol)
dy8f1 = M(1)*yp(1) + M(2)*yp(2) + M(3)*yp(3) + M(4)*yp(4) + M(5)*yp(5);
dy8f2 = 28*yp(1) + 28*yp(2) + 28*yp(3) + 28*yp(4) + 28*yp(5);
dy(8) = dy8f1 - dy8f2 - yp(8);

% -------------------------------------------------------------------------
% Pressure: dy(9) ---------------------------------------------------------
% -------------------------------------------------------------------------

% Constants
mu = 0.048; % lb/(ft-h)
mus = mu/3600; %lb(ft-s)
eB = 0.5;
dpft = 0.25/12; % ft
Rg = 19.314; % psia ft3/lbmol/K

% Additional Functions
alpha = 1 + (2*dpft)/(3*(1-eB)*dt); %unitless
u = ( (359)*(14.7)*W*y(6) ) / ( (273)*(3600)*(NT)*Ac*y(8)*y(9) );  %ft/s

% dy(9) = dPT (Psia)
dy9f1 = ( -alpha^2 ) / ( 32.2*144 );
dy9f2 = ( (150*mus*u)*((1-eB)^2) ) / ( (dpft^2)*(eB^3) );
dy9f3 = ( 1.75*y(8)*(u^2)*(1-eB)*y(9) ) / ( alpha*Rg*y(6)*dpft*((eB)^3) );
dy(9) = ( dy9f1 * (dy9f2 + dy9f3) ) - yp(9);

function dy = systemeqs2(z,y,yp)

% -------------------------------------------------------------------------
% Countercurrent System (mode = -1) ---------------------------------------
% -------------------------------------------------------------------------

% Define f(z,y,yp)=0 vector, implicit diff eq system
dy = ones(9,1);

% -------------------------------------------------------------------------
% Rates (lbmol / lbcat h) -------------------------------------------------
% -------------------------------------------------------------------------

% RT and Partial Pressure Expressions
RT = ( (1/y(6))-(1/373) )/1.987;
PDT = y(1)*y(9); %psia
PVA = y(2)*y(9); %psia
PO2 = y(3)*y(9); %psia

% Reaction Constants (k1-k12)
k = [ 1.771*10^-3; 23295.0; 0.5; 1.0; 0.8184; 0.5; 0.2314;
    1.0; 1.25; 2.0; 2.795*10^-4; 33000.0; 0.5; 2.0];

% Reaction 1 Rate Expression
r1top = (k(1)*exp(-k(2)*RT)*(PO2^k(3))*(PDT^k(4)));
r1bom =(1 + k(5)*(PO2^k(6))+ k(7)*(PDT^k(8))+ k(9)*(PVA^k(10)));
ra = r1top / r1bom; % lbmol/(lbcat-h)

% Reaction 2 Rate Expresion
r2top = (k(11)*exp(-k(12)*RT)*(PO2^k(13))*(PVA^k(14)));
r2bom =(1 + k(5)*(PO2^k(6))+ k(7)*(PDT^k(8))+ k(9)*(PVA^k(10)));
rb = r2top / r2bom; % lbmol/(lbcat-h)

R1 = -ra;
R2 = ra - rb;
R3 = -0.5*ra - 2.5*rb;
R4 = 2*rb;
R5 = ra + 2*rb;

% -------------------------------------------------------------------------
% Mass Balances: dy(1)-dy(5)  ---------------------------------------------
% -------------------------------------------------------------------------

% Constants
NT = 2500; % number tubes
rhoB = 100; %lb/ft3
dt = 1/12; % ft
rt = dt/2; % ft
Ac = pi*(rt^2); %ft^2
W = 100000; % lb/h

% dy(1) = dX1, X1-X5 unitless
dy1f1 = ((NT*rhoB*Ac*y(8))/W)*(R1);
dy1f2 = (y(1)/y(8))*yp(8);
dy(1) = dy1f1 + dy1f2 - yp(1);

% dy(2) = dX2, X1-X5 unitless
dy2f1 = ((NT*rhoB*Ac*y(8))/W)*(R2);
dy2f2 = (y(2)/y(8))*yp(8);
dy(2) = dy2f1 + dy2f2 - yp(2);

% dy(3) = dX3, X1-X5 unitless
dy3f1 = ((NT*rhoB*Ac*y(8))/W)*(R3);
dy3f2 = (y(3)/y(8))*yp(8);
dy(3) = dy3f1 + dy3f2 - yp(3);

% dy(4) = dX4, X1-X5 unitless
dy4f1 = ((NT*rhoB*Ac*y(8))/W)*(R4);
dy4f2 = (y(4)/y(8))*yp(8);
dy(4) = dy4f1 + dy4f2 - yp(4);

% dy(5) = dX5, X1-X5 unitless
dy5f1 = ((NT*rhoB*Ac*y(8))/W)*(R5);
dy5f2 = (y(5)/y(8))*yp(8);
dy(5) = dy5f1 + dy5f2 - yp(5);

% -------------------------------------------------------------------------
% Energy Balances: dy(6)-dy(7) --------------------------------------------
% -------------------------------------------------------------------------

% Constants
mode = -1; % type of flow
U = 120.5; % BTU/(h-ft^2-C)
FcCpc = 10^6; % BTU/(h C)
Hr1 = -74.32*1000; % BTU/lbmol
Hr2 = -474.57*1000; % BTU/lbmol
Cp = 0.50; % BTU/(lb C)

% dy(6) = dTR, TR (K)
dy6f1 = (-1*NT)*Ac*(Hr1*ra + Hr2*rb)*rhoB;
dy6f2 = pi*NT*dt*U*( y(6)-y(7) );
dy6f3 = W*Cp;
dy(6) = ((dy6f1 - dy6f2)/dy6f3) - yp(6);

% dy(7) = dTC, TC (K)
dy(7) = ((mode)*pi*NT*dt*U*(y(6)-y(7))) / ((FcCpc)) - yp(7);

% -------------------------------------------------------------------------
% Molecular Weight: dy(8) -------------------------------------------------
% -------------------------------------------------------------------------

% Molecular Weights
M = [46; 44; 32; 44; 18]; % lb/lbmol

% dy(8) = dMW (lb/lbmol)
dy8f1 = M(1)*yp(1) + M(2)*yp(2) + M(3)*yp(3) + M(4)*yp(4) + M(5)*yp(5);
dy8f2 = 28*yp(1) + 28*yp(2) + 28*yp(3) + 28*yp(4) + 28*yp(5);
dy(8) = dy8f1 - dy8f2 - yp(8);

% -------------------------------------------------------------------------
% Pressure: dy(9) ---------------------------------------------------------
% -------------------------------------------------------------------------

% Constants
mu = 0.048; % lb/(ft-h)
mus = mu/3600; %lb(ft-s)
eB = 0.5;
dpft = 0.25/12; % ft
Rg = 19.314; % psia ft3/lbmol/K

% Additional Functions
alpha = 1 + (2*dpft)/(3*(1-eB)*dt); %unitless
u = ( (359)*(14.7)*W*y(6) ) / ( (273)*(3600)*(NT)*Ac*y(8)*y(9) );  %ft/s

% dy(9) = dPT (Psia)
dy9f1 = ( -alpha^2 ) / ( 32.2*144 );
dy9f2 = ( (150*mus*u)*((1-eB)^2) ) / ( (dpft^2)*(eB^3) );
dy9f3 = ( 1.75*y(8)*(u^2)*(1-eB)*y(9) ) / ( alpha*Rg*y(6)*dpft*((eB)^3) );
dy(9) = ( dy9f1 * (dy9f2 + dy9f3) ) - yp(9);
