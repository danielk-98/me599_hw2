syms theta_1
syms theta_2  
syms theta_3
syms theta_1dot
syms theta_2dot
syms theta_3dot




L1 = 0.4;
L2 = 1.0;
LC2 = 0.5;
LC1 = 0.2;
LC3 = 0.5;
rC3 = 0.1;

I = [1, 0, 0; 0, 1, 0; 0, 0, 1]; 
zero = [0;0;0];

T01 = [cos(theta_1), -sin(theta_1), 0, 0;
    sin(theta_1), cos(theta_1), 0, 0;
    0, 0, 1, 0;
    0,0,0,1 ]

T1C1(1:3, 1:3) = I; 
T1C1(1:3, 4) = [0;0;LC1];
T1C1(4, 1:4) = [0,0,0,1];

T1C1

T12(1:3, 1:3) = I;
T12(1:3, 4) = [0; 0; L1];
T12(4, 1:4) = [0,0,0,1];

T12 = T12 * [1,0,0,0;0,0,-1,0;0,1,0,0;0,0,0,1]; 
T12 = T12 * [cos(theta_2), -sin(theta_2), 0, 0; 
            sin(theta_2), cos(theta_2), 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1;];
T12

T2C2(1:3, 1:3) = I;
T2C2(1:3, 4) = [0;LC2;0];
T2C2(4, 1:4) = [0,0,0,1]; 

T2C2

T23(1:3, 1:3) = I;
T23(1:3, 4) = [0; L2; 0];
T23(4, 1:4) = [0,0,0,1];
T23 = T23 * [cos(theta_3), -sin(theta_3), 0, 0;
            sin(theta_3), cos(theta_3), 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1;];
T23


T3C3(1:3, 1:3) = I;
T3C3(1:3, 4) = [0; LC3; rC3];
T3C3(4, 1:4) = [0,0,0,1];

T3C3

R01 = T01(1:3, 1:3);
T0C1 = T01*T1C1;
oc1 = T0C1(1:3, 4); 
o1 = T01(1:3,4);

z = [0;0;1];

JC1 = [cross([R01 * z], oc1 - o1), zero, zero ; R01*z, zero, zero]

T02 = T01*T12
R02 = T02(1:3, 1:3); 
T0C2 = T01*T12*T2C2; 
oc2 = T0C2(1:3, 4);
o2 = T02(1:3, 4);

JC2 = [cross([R01 * z], oc2 - o1), cross([R02 * z], oc2 - o2), zero; R01 * z, R02*z, zero]

T03 = T01*T12*T23; 
R03 = T03(1:3, 1:3); 
T0C3 = T03*T3C3; 
oc3 = T0C3(1:3, 4); 
o3 = T03(1:3, 4);

JC3 = [cross([R01 * z], oc3 - o1), cross([R02 * z], oc3 - o2), cross([R03*z], oc3 - o3); 
    R01*z, R02*z, R03*z]

qdot = [theta_1dot;theta_2dot;theta_3dot];

I1 = [1,0,0;0,0.083,0;0,0,1;];
m1 = 1;
KE1 = 1/2*transpose(qdot)*(transpose(JC1(4:6,:))*(R01*I1*transpose(R01))*JC1(4:6,:) + m1*transpose(JC1(1:3,:))*JC1(1:3,:))*qdot


I2 = [1,0,0;0,0.083,0;0,0,1;];
m2 = 1;
KE2 = 1/2*transpose(qdot)*(transpose(JC2(4:6,:))*(R02*I2*transpose(R02))*JC2(4:6,:) + m2*transpose(JC2(1:3,:))*JC2(1:3,:))*qdot

I3 = [1,0,0;0,0.33,0;0,0,1;];
m3 = 1;
KE3 = 1/2*transpose(qdot)*(transpose(JC3(4:6,:))*(R03*I3*transpose(R03))*JC3(4:6,:) + m3*transpose(JC3(1:3,:))*JC3(1:3,:))*qdot

KE = KE1 + KE2 + KE3;

g = 9.81;

PE1 = m1*g*T0C1(3,4);
PE2 = m2*g*T0C2(3,4);
PE3 = m3*g*T0C3(3,4); 


PE = PE1 + PE2 + PE3;

L = KE - PE

dLdqdot1 = diff(L, theta_1dot);
dLdqdot2 = diff(L, theta_2dot); 
dLdqdot3 = diff(L, theta_3dot);
dLdq1 = diff(L, theta_1); 
dLdq2 = diff(L, theta_2); 
dLdq3 = diff(L, theta_3); 
% pretty(dLdQ)


