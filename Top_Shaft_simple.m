function [ Defl theta ] = Top_Shaft_simple( D_shaft, Bearing_center_distance )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


d_pulley=34; % shaft diameter for the pulley 
d_bearing=35; % bearing bore diameter in mm
d_shaft=41; % min shaft shoulder diameter in mm
bear_thick=17; % bearing thickness in mm
pulley_thick=40; %width of pulley in mm
spl= 300; %spool length in meters
g=9.81; %gravity in L/s^2
F_spl=50*g; % weight of spool and wire just before the winding is complete
F_pul=500;  % resultant pulley tension

Gear_ratio=0.6;
T1=30/pi; % Nm torque of the motor

T2=0.94*T1/Gear_ratio;

T_wire=T2/0.12; % 12 cm wire winding radius when spool is almost full


w= 76.5e-6;  % carbon steel(AISI 1050) unit weight in N/mm^3 (table A-5)
E=207e9; % carbon steel(AISI 1050) modulus of elasticity in Pa (table A-5)

if D_shaft>=41
    d_shaft=D_shaft;
end

if Bearing_center_distance>bear_thick
    S=Bearing_center_distance;
else
    S=bear_thick; % minimum bearing center distance
end


%I_inv=Inertia_inv(d_shaft,S);


W1= w*pi*d_bearing^2;  % weight of shaft segment near the spool and ...
                            %pulley in N/mm
                            
W2= w*pi*d_shaft^2;  % weight of the shaft's middle segment 

W3=W1;

W4= w*pi*d_pulley^2; 



% H=bear_thick/2; %8.5
% I=spl/2; %150
% J=spl+H; %308.5
% K=spl+bear_thick; %317 first shoulder position
% L=spl+S; %300+S
% 
% length= K+S+40;

% syms pos C0 C1

% V=@(x) -W1*x*step2(x)-F_spl*step2(x-I)+Ra_y*step2(x-J)-...
%         (W2-W1)*(x-K)*step2(x-K)-(W1-W2)*(x-L)*step2(x-L)+...
%          Rb_y*step2(x-(L+H))-F_pul*step2(x-(K+S+10)); %shear in Newtons
% 
%     
% M= matlabFunction(int(V(pos)));

% V=@(x) -W1*x*step2(x)-F_spl*step2(x-I)+Ra_y*step2(x-J)-...
%         (W2-W1)*(x-K)*step2(x-K)-(W1-W2)*(x-L)*step2(x-L)+...
%          Rb_y*step2(x-(L+H))-F_pul*step2(x-(K+S+10)); %shear in Newtons

    
%M= matlabFunction(int(V(pos)));
   
  
%theta_prime= @(x)(1/E)*M(x)*I_inv(x);

% all units in mm

syms pos


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
                            dfspl,K,L,N,O,T_wire,F_spl,F_pul);


V=@(x) -W1*x*step2(x)-F_spl*step2(x-dfspl)+Ra_y*step2(x-db1)-...
        (W2-W1)*(x-K)*step2(x-K)-(W3-W2)*(x-L)*step2(x-L)+...
         Rb_y*step2(x-db2)-(W4-W3)*(x-N)*step2(x-N)-...
         F_pul*step2(x-dfpul); %shear in Newtons
     
M=matlabFunction(int(V(pos)));


I1=(pi/64)*(0.035)^4;
I2=(pi/64)*(0.045)^4;
I4=(pi/64)*(0.034)^4;


theta1=matlabFunction(int(M(pos))/E*I1);
theta2=matlabFunction(int(M(pos))/E*I2);
theta2= @(x) theta2(x)+theta1(K)-theta2(K);
theta3=theta1;
theta3= @(x) theta3(x)+theta2(L)-theta3(L);
theta4=matlabFunction(int(M(pos))/E*I4);
theta4= @(x) theta4(x)+theta3(N)-theta4(N);

theta=@(x)theta1(x)*(1-step2(x-K))+theta2(x)*(step2(x-K)-step2(x-L))+...
    theta3(x)*(step2(x-L)-step2(x-N))+theta4(x)*step2(x-N);


def1=matlabFunction(int(int(M(pos)))/E*I1);
def2=matlabFunction(int(int(M(pos)))/E*I2);
def2= @(x) def2(x)+(theta1(K)-theta2(K))*x+def1(K)-def2(K);
def3=def1;
def3=@(x) def3(x)+(theta2(L)-theta3(L))*x+def2(L)-def3(L);
def4=matlabFunction(int(int(M(pos)))/E*I4);
def4=@(x) def4(x)+(theta3(N)-theta4(N))*x+def3(N)-def4(N);

Defl=@(x)def1(x)*(1-step2(x-K))+def2(x)*(step2(x-K)-step2(x-L))+...
    def3(x)*(step2(x-L)-step2(x-N))+def4(x)*step2(x-N);


% theta_var=matlabFunction(int(M(pos))*I_inv(pos));
% 
% 
% Defl_var= matlabFunction(int(int(M(pos)))*I_inv(pos));
% 
% 
% def_b1=Defl_var(J)/(E*10^6)+C1*J+C0;
% def_b2=Defl_var(J+S)/(E*10^6)+C1*(J+S)+C0;
% 
% Coeff=solve([def_b1==0,def_b2==0]);
% 
% C0_const=Coeff.C0;
% C1_const=Coeff.C1;
% 
% 
% theta=@(x) double((theta_var(x)/(E*10^6))+C1_const);%in rad
% 
% Defl=@(x) double((Defl_var(x)/(E*10^6))+C1_const*x+C0_const); %deflection in mm


    
end

