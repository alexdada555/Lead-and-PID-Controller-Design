clear;
clc;
format long;

%% SYSTEM UNDER CONSIDERATION

G = tf(822,[1 14 63 0]);
% Second order system with integrator
figure(1)
rlocus(G);
%% Requirements Frequency Domain Calculations
% Percentage overshoot < 15%
% Rise Time < 100ms
% Settling Time < 500ms within 5%

A = log(15);
Rzeta = A/(sqrt(pi^2 + A^2)); % desired damping ratio
CPm = 100*Rzeta; % desired compensated phase magring in degrees

% tan(Cpm) = 8/(Ts*Wc)
%increasing threshold margin
% At Rwc of 7.36, Phase angle = -175
% However i will aim at a number < Ts of 200ms
% to ensure target is met
Rwc = (8)/(0.20*tan(deg2rad(CPm))); 

% At Rwc of 18.4056, Phase angle = -227
% additional Phase angle is:
AdPa = CPm + 227 - 180;

%% LEAD COMPENSATOR POLE AND ZERO CALCULATIONS

CPm1 = deg2rad(AdPa+5); % added 5 degrees
% sgrt(Z*P) = Rwc
alpha = ((1 + sin(CPm1))/(1 - sin(CPm1)));

zs = Rwc/sqrt(alpha);
ps = alpha*zs;

% lead compensator is
Gcl = tf([1/zs 1],[1/ps 1]);

%% PLOT OF SINGLE LEAD COMPENSATED SYSTEM

Xc = ((G*Gcl)/(1+(G*Gcl)));

% A single lead compensator proves to be insufficient to
% meet the system requirements only providing an additional
% 30 degrees of phase margin

%% PLOT OF DOUBLE LEAD COMPENSATED SYSTEM

Xdoub = ((G*Gcl*Gcl)/(1+(G*Gcl*Gcl)));
% A double lead compensator is used

figure(2);
step(Xc);
stepdataA = stepinfo(Xc,'SettlingTimeThreshold',0.05);
hold all;
step(Xdoub);
stepdataB = stepinfo(Xdoub,'SettlingTimeThreshold',0.05);
legend('Single Lead','Double Lead');
grid on;

figure(3);
margin(G);
hold all;
margin(Xc);
hold all;
margin(Xdoub);
legend('Original','Single Lead','Double Lead');
grid on;

%% DISCRETISATION OF SYSTEM AND ZtoW-DOMAIN 

% New Desired cross over frequency in  W-domain
% With New Desired Ts at 300ms instead of 200ms
% used in  S-domain. obtained through trial and error
Rwcw = (8)/(0.30*tan(deg2rad(CPm)));
% The desired system bandwidth
ft = ((4/0.30)*(1/Rzeta)*sqrt(((1-(2*(Rzeta^2)))+sqrt((4*(Rzeta^4))-(4*(Rzeta^2))+2))));
T = 1/(20*ft);
% In order to determine the sampling time the system
% bandwidth is calculated and 20 X bandwidth is used
% to avoid interference and aliasing.
[ZZnum,ZZdnum] = c2dm(822,[1 14 63 0],T);
% The Z-transform of G is taken
minZZnum = [-1*ZZnum(1) ZZnum(2) -1*ZZnum(3) ZZnum(4)];
minZZdnum = [-1*ZZdnum(1) ZZdnum(2) -1*ZZdnum(3) ZZdnum(4)];

[YYnum,YYdnum] = bilinear(minZZnum,minZZdnum,0.5);
% Due to the nature of the bilinear transformation
% -Z is used instead of Z
factor = 1/((-T/2)^3);
WWnum = YYnum.*[(-T/2)^3 (-T/2)^2 -T/2 1]*(factor);
WWdnum = YYdnum.*[(-T/2)^3 (-T/2)^2 -T/2 1]*(factor);
Gw = tf(WWnum,WWdnum);
% In order to translate the bilinear transform into
% the W-domain  each coefficient 1 to 1 multplied by 
% increaing powers of -T/3

%% W-DOMAIN LEAD COMPENSATOR POLE AND ZERO CALCULATIONS
% At Rwcw of 12.2704, Phase angle = 153
% additional Phase angle is:
AdPaw = CPm + 153 - 180;

CPm1w = deg2rad(AdPaw+5); % added 5 degrees
% sgrt(Z*P) = Rwc
alphaw = ((1 + sin(CPm1w))/(1 - sin(CPm1w)));

zw = Rwc/sqrt(alphaw);
pw = alphaw*zw;

% lead compensator is
Dw = tf([1/zw 1],[1/pw 1]);

%% PLOT OF W-DOMAIN SINGLE LEAD COMPENSATED SYSTEM

Xcw = ((Gw*Dw*Dw*Dw)/(1+(Gw*Dw*Dw*Dw)));

figure(4);
margin(Gw);
hold all
margin(Xcw);
legend('W-Domain','W-Triple Lead Compensated');
grid on;

figure(5);
step(Xdoub);
hold all;
step(Xcw);
stepdataC = stepinfo(Xcw,'SettlingTimeThreshold',0.05);
legend('S-domain','W-domain');
grid on;

%% CONVERSION FROM WtoZ-DOMAIN
z = tf('z',T);

Wzero = [1/zw 1];
Wpole = [1/pw 1];

w = (((z-1)/(z+1))*(2/T));
%z = (w*T/2)*(1-w)/(w+1)
Dznum = ((w/zw)+1);
Dzdnum = ((w/pw)+1);
%Dznum = (((2*z)-2+(z*zw*T)+(zw*T))/((z*zw*T)+(zw*T)));
%Dzdnum = (((2*z)-2+(z*pw*T)+(pw*T))/((z*pw*T)+(pw*T)));

Dz = ((Dznum)/(Dzdnum));
%(((w/zw)+1)/((w/zw)+1));

Gz = tf(ZZnum,ZZdnum,T);
XCz = feedback(Gz,Dz*Dz);
%XCz = ((Gz*Dz)/(1+(Gz*Dz)));

figure(6);
step(XCz);
