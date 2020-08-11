A = -0.1976068851;
K = 198.11;
B = 0.205112;
C = 0.557971;
D = 1;

Gp = tf([A*K K],[B C D]);

zeta = 0.4559;
Wn = 2.3532;
P = 10*zeta*Wn;

syms Kp Ki Kd
eq1 = C/(B - A*K*Kd) + (K*Kd)/(B - A*K*Kd) - (A*K*Kp)/(B - A*K*Kd) == P + 2*zeta*Wn;
eq2 = D/(B - A*K*Kd) + (K*Kp)/(B - A*K*Kd) - (A*K*Kp)/(B - A*K*Kd) == 2*Wn*P + Wn^2;
eq3 = (A*K*Ki)/(B - A*K*Kd) == P*(Wn^2);

S = solve(eq1,eq2,eq3);

X = double(vpa(S.Kp))
Y = double(vpa(S.Ki))
Z = double(vpa(S.Kd))

Gc = tf([Z X Y],[1 0]);
%Gc = tf([0.005 0.002 0.02],[1 0]);
Gc = tf([0.0035 0.01 0.02],[1 0]);
%Gc = tf([1 5],[1 10]);

figure(1);
%step(Gp);
%hold on;
GK = feedback(Gp,Gc);
step(GK);

figure(2);
Gc = tf([0.0005 0.0005 0.0208],[(0.0035/ 0.009)/20 1]);
GK = feedback(Gp,Gc);
step(GK);
