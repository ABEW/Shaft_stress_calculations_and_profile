function [ theta deflection M V] = Top_shaft_final( D_shaft, Bearing_center_distance )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

d_pulley=34; % shaft diameter for the pulley 
d_bearing=35; % bearing bore diameter in mm
d_shaft=41; % min shaft shoulder diameter in mm
bear_thick=17; % bearing thickness in mm
pulley_thick=40; %width of pulley in mm
spl= 300; %spool length in meters
g=9.81; %gravity in m/s^2
F_spl=50*g; % weight of spool and wire just before the winding is complete
F_pul=500;  % resultant pulley tension

Gear_ratio=0.6;
T1=30/pi; % Nm torque of the motor

T2=0.94*T1/Gear_ratio;

T_wire=T2/0.12; % 12 cm wire winding radius when spool is almost full


w= 76.5e-6;  % carbon steel(AISI 1050) unit weight in N/mm^3 (table A-5)
E=207e9; % carbon steel(AISI 1050) modulus of elasticity in Pa (table A-5)

if D_shaft>=41
    d_shaft=41;
    'shaft/shoulder diameter set to default value of 41 mm'
end

if Bearing_center_distance>bear_thick
    S=Bearing_center_distance;
else
    S=bear_thick; % minimum bearing center distance
    'bearing distance set to default value of 17mm center to center'
end


W1= w*pi*d_bearing^2;  % weight of shaft segment near the spool in N/mm
                            
W2= w*pi*d_shaft^2;  % weight of the shaft's middle segment 

W3=W1; %weight of shaft segment near bearing 2  in N/mm

W4= w*pi*d_pulley^2;  %weight of shaft segment near pulley in N/mm


% all units in mm

syms pos C0 C1


clearance=3; %clearance between the bearing and spool and bearing-pulley
db1=300+clearance+bear_thick/2;  %distance to the first bearing from left end of shaft
db2=db1+S; %distance to the second bearing from left end of shaft
dfspl=150;  % distance to the center of spool
K=db1+bear_thick/2; %distance to the first shoulder
L=K+S-bear_thick; %distance to the second shoulder
N=L+bear_thick; %distance to the third shoulder
O=N+clearance+pulley_thick; % total length of the shaft
dfpul=O-30; %distance to the center of belt (not pulley)

[Ra_y,Rb_y,Ra_z,Rb_z]=Reaction_forces2(W1,W2,W3,W4,db1,db2,dfpul,...
                            dfspl,K,L,N,O,T_wire,F_spl,F_pul)


V=@(x) -W1*x*step2(x)-F_spl*step2(x-dfspl)+Ra_y*step2(x-db1)-...
        (W2-W1)*(x-K)*step2(x-K)-(W3-W2)*(x-L)*step2(x-L)+...
         Rb_y*step2(x-db2)-(W4-W3)*(x-N)*step2(x-N)-...
         F_pul*step2(x-dfpul); %shear in Newtons
     
M=matlabFunction(int(V(pos)));


I1=(pi/64)*(0.035)^4;
I2=(pi/64)*(0.045)^4;
I3=I1;
I4=(pi/64)*(0.034)^4;

% used to test if the M_over_I function below is accurate

% I_x= @(x) I1*step2(x-0)+(I2-I1)*step2(x-K)+... 
%     (I3-I2)*step2(x-L)+(I4-I3)*step2(x-N);

%test_ratio=@(x) M(x)/I_x(x);
%

g= @(x) -0.5*W1*x^2-F_spl*(x-dfspl)+Ra_y*(x-db1)-0.5*(W2-W1)*(x-K)^2;
    % just a function I used to simplify the below function


M_over_I= @(x) (-0.5*W1*x^2*step2(x)-F_spl*(x-dfspl)*step2(x-dfspl)+Ra_y*(x-db1)*step2(x-db1))/I1+...
    ((1/I2)-(1/I1))*step2(x-K)*(g(x)+0.5*(W2-W1)*(x-K)^2)-((W2-W1)*(x-K)^2*step2(x-K))/(2*I2)+...
    ((1/I3)-(1/I2))*step2(x-L)*g(x)+((-0.5*(W3-W2)*(x-L)^2*step2(x-L))+Rb_y*(x-db2)*step2(x-db2))/I3+...
    ((1/I4)-(1/I3))*step2(x-N)*(g(x)-0.5*(W3-W2)*(x-L)^2+Rb_y*(x-db2))+...
    (-1/I4)*(0.5*(W4-W3)*(x-N)^2*step2(x-N)+F_pul*(x-dfpul)*step2(x-dfpul)); %in N.mm/m^4

theta_prime=@(x) (1/E)*M_over_I(x); % d^y/dx^2 in mm/m^2

theta1=matlabFunction(int(theta_prime(pos))); % dy/dx in mm^2/m^2 

theta_var=@(x) (1/10^6)*theta1(x); % dy/dx in radians w/o constants

deflection_var=matlabFunction(int(theta_var(pos))); % in mm w/o constants

% The steps below are used to determine the integration coefficients using
% the restriction that deflection at the bearings must be zero

def_b1=deflection_var(db1)+C1*db1+C0;
def_b2=deflection_var(db2)+C1*db2+C0;

Coeff=solve([def_b1==0,def_b2==0]);
 
C0_const=double(Coeff.C0);
C1_const=double(Coeff.C1);

theta= @(x) theta_var(x)+C1_const; % slope in radians
deflection=@(x) deflection_var(x)+C1_const*x+C0_const; % deflection in mm
    
    
end

