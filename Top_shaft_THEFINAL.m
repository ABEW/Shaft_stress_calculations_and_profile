function [ theta_Y theta_Z deflection_Y deflection_Z M_Y M_Z] = Top_shaft_THEFINAL( D_shaft, Bearing_center_distance )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all;

d_pulley=34; % shaft diameter for the pulley
d_bearing=35; % bearing bore diameter in mm

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

if D_shaft<41
    d_shaft=D_shaft;
else
    d_shaft=41; % min shaft shoulder diameter in mm
    'shaft/shoulder diameter set to default value of 41 mm'
end

if Bearing_center_distance>bear_thick
    S=Bearing_center_distance+bear_thick;
else
    S=2*bear_thick; % minimum bearing center distance
    'bearing inner edge distance set to default value of 17mm'
end


W1= 0.25*w*pi*d_bearing^2;  % weight of shaft segment near the spool in N/mm

W2= 0.25*w*pi*d_shaft^2;  % weight of the shaft's middle segment

W3=W1; %weight of shaft segment near bearing 2  in N/mm

W4= 0.25*w*pi*d_pulley^2;  %weight of shaft segment near pulley in N/mm


% all units in mm

syms pos C0_Y C1_Y C0_Z C1_Z


clearance=3; %clearance between the bearing and spool and bearing-pulley
db1=300+bear_thick/2;  %distance to the first bearing from left end of shaft
db2=db1+S; %distance to the second bearing from left end of shaft
dfspl=150;  % distance to the center of spool
K=db1+bear_thick/2; %distance to the first shoulder
L=K+S-bear_thick; %distance to the second shoulder
N=L+bear_thick; %distance to the third shoulder
 ' total length of the shaft'
O=N+clearance+pulley_thick
dfpul=O-30; %distance to the center of belt (not pulley)

[Ra_y,Rb_y,Ra_z,Rb_z]=Reaction_forces2(W1,W2,W3,W4,db1,db2,dfpul,...
    dfspl,K,L,N,O,T_wire,F_spl,F_pul);

Torsion=@(x) T2*(step2(x-10)-step2(x-dfpul));

Figuretastic(Torsion,[0,O+10],'Torque diagram','Distance from edge near the spool(x) [mm]','|Torque| [Nm]')


V_Y=@(x) -W1*x*step2(x)-F_spl*step2(x-dfspl)+Ra_y*step2(x-db1)-...
    (W2-W1)*(x-K)*step2(x-K)-(W3-W2)*(x-L)*step2(x-L)+...
    Rb_y*step2(x-db2)-(W4-W3)*(x-N)*step2(x-N)-...
    F_pul*step2(x-dfpul)+W4*(x-O)*step2(x-O); %Y direction shear in Newtons

Figuretastic(V_Y,[0,O+10],'Shear force (Y-direction) diagram','Distance from edge near the spool(x) [mm]','Shear Force (V_y) [N]')


M_Z=matlabFunction(int(V_Y(pos))); % Z-moment from V_Y in N mm
Figuretastic(@(x)M_Z(x)/10^3,[0,O+10],'Moment (Z-direction) diagram','Distance from edge near the spool(x) [mm]','Moment (M_z) [Nm]')


V_Z=@(x) T_wire*step2(x)-Ra_z*step2(x-db1)-Rb_z*step2(x-db2); %Z-direction shear in Newtons
Figuretastic(V_Z,[0,O+10],'Shear force (Z-direction) diagram','Distance from edge near the spool(x) [mm]','Shear Force (V_z) [N]')


M_Y=matlabFunction(int(-1*V_Z(pos))); % Y-moment from V_Y in N mm
Figuretastic(@(x)M_Y(x)/10^3,[0,O+10],'Moment (Y-direction) diagram','Distance from edge near the spool(x) [mm]','Moment (M_y) [Nm]')



I1=(pi/64)*(0.035)^4;
I2=(pi/64)*(0.045)^4;
I3=I1;
I4=(pi/64)*(0.034)^4;

% used to test if the M_over_I_Y function below is accurate

% I_x= @(x) I1*step2(x-0)+(I2-I1)*step2(x-K)+...
%     (I3-I2)*step2(x-L)+(I4-I3)*step2(x-N);

%test_ratio=@(x) M(x)/I_x(x);
%

G= @(x) -0.5*W1*x^2-F_spl*(x-dfspl)+Ra_y*(x-db1)-0.5*(W2-W1)*(x-K)^2;
% just a function I used to simplify the M/I_Y function

H= @(x) -T_wire*x+Ra_z*(x-db1);


M_over_I_Y= @(x) (-0.5*W1*x^2*step2(x)-F_spl*(x-dfspl)*step2(x-dfspl)+Ra_y*(x-db1)*step2(x-db1))/I1+...
    ((1/I2)-(1/I1))*step2(x-K)*(G(x)+0.5*(W2-W1)*(x-K)^2)-((W2-W1)*(x-K)^2*step2(x-K))/(2*I2)+...
    ((1/I3)-(1/I2))*step2(x-L)*G(x)+((-0.5*(W3-W2)*(x-L)^2*step2(x-L))+Rb_y*(x-db2)*step2(x-db2))/I3+...
    ((1/I4)-(1/I3))*step2(x-N)*(G(x)-0.5*(W3-W2)*(x-L)^2+Rb_y*(x-db2))+...
    (-1/I4)*(0.5*(W4-W3)*(x-N)^2*step2(x-N)+F_pul*(x-dfpul)*step2(x-dfpul)-0.5*W4*(x-O)^2*step2(x-O)); %in N.mm/m^4


M_over_I_Z= @(x) (-T_wire*x*step2(x)+Ra_z*(x-db1)*step2(x-db1)+Rb_z*(x-db2)*step2(x-db2))/I1+...
    (((1/I2)-(1/I1))*step2(x-K)+((1/I3)-(1/I2))*step2(x-L))*H(x)+...
    ((1/I3)-(1/I1))*Rb_z*(x-db2)*step2(x-db2)+...
    ((1/I4)-(1/I3))*step2(x-N)*(H(x)+Rb_z*(x-db2));


theta_prime_Y=@(x) (1/E)*M_over_I_Y(x); % d^y/dx^2 in mm/m^2

theta_prime_Z=@(x) (1/E)*M_over_I_Z(x); % d^z/dx^2 in mm/m^2


theta1_Y=matlabFunction(int(theta_prime_Y(pos))); % dy/dx in mm^2/m^2

theta1_Z=matlabFunction(int(theta_prime_Z(pos))); % dz/dx in mm^2/m^2


theta_var_Y=@(x) (1/10^6)*theta1_Y(x); % dy/dx in radians w/o constants

theta_var_Z=@(x) (1/10^6)*theta1_Z(x); % dz/dx in radians w/o constants


deflection_var_Y=matlabFunction(int(theta_var_Y(pos))); % in mm w/o constants

deflection_var_Z=matlabFunction(int(theta_var_Z(pos))); % in mm w/o constants

% The steps below are used to determine the integration coefficients using
% the restriction that deflection at the bearings must be zero

def_b1_Y=deflection_var_Y(db1)+C1_Y*db1+C0_Y;
def_b2_Y=deflection_var_Y(db2)+C1_Y*db2+C0_Y;

def_b1_Z=deflection_var_Z(db1)+C1_Z*db1+C0_Z;
def_b2_Z=deflection_var_Z(db2)+C1_Z*db2+C0_Z;

Coeff_Y=solve([def_b1_Y==0,def_b2_Y==0]);
Coeff_Z=solve([def_b1_Z==0,def_b2_Z==0]);

C0_Y_const=double(Coeff_Y.C0_Y);
C1_Y_const=double(Coeff_Y.C1_Y);

C0_Z_const=double(Coeff_Z.C0_Z);
C1_Z_const=double(Coeff_Z.C1_Z);

theta_Y= @(x) theta_var_Y(x)+C1_Y_const; % slope in radians
deflection_Y=@(x) deflection_var_Y(x)+C1_Y_const*x+C0_Y_const; % deflection in mm

Figuretastic(theta_Y,[0,O+10],'Shaft Slope (Y-direction)','Distance from edge near the spool(x) [mm]','Slope (\theta_y) [rad]')
Figuretastic(deflection_Y,[0,O+10],'Shaft Deflection (Y-direction)','Distance from edge near the spool(x) [mm]','Deflection(y) [mm]')


theta_Z= @(x) theta_var_Z(x)+C1_Z_const; % slope in radians
deflection_Z=@(x) deflection_var_Z(x)+C1_Z_const*x+C0_Z_const; % deflection in mm

Figuretastic(theta_Z,[0,O+10],'Shaft Slope (Z-direction)','Distance from edge near the spool(x) [mm]','Slope (\theta_x) [rad]')
Figuretastic(deflection_Z,[0,O+10],'Shaft Deflection (Z-direction)','Distance from edge near the spool(x) [mm]','Deflection(x) [mm]')




fprintf('\nDEFLECTION AT REQUIRED POINTS \n')

def_spool=norm([deflection_Z(dfspl),deflection_Y(dfspl)]);
def_pulley=norm([deflection_Z(dfpul),deflection_Y(dfpul)]);

fprintf('\nDeflection at spool is: %f mm \n',def_spool)

fprintf('Deflection at pulley is: %f mm \n',def_pulley)



fprintf('\nSLOPE AT REQUIRED POINTS \n')

slope_spool=norm([theta_Z(dfspl),theta_Y(dfspl)]);
slope_pulley=norm([theta_Z(dfpul),theta_Y(dfpul)]);
slope_bearing1=norm([theta_Z(db1),theta_Y(db1)]);
slope_bearing2=norm([theta_Z(db2),theta_Y(db2)]);

fprintf('\nSlope at spool is: %f radians \n',slope_spool)

fprintf('Slope at pulley is: %f radians \n',slope_pulley)

fprintf('Slope at left bearing is: %f radians \n',slope_bearing1)

fprintf('Slope at right bearing is: %f radians \n',slope_bearing2)
fprintf('\n')


end

