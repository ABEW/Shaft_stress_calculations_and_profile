function [ Defl2 Slope2 I_x ] = Top_Shaft( D_shaft, Bearing_center_distance )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Variables
d_bearing=35; % bearing bore diameter in mm
d_shaft=41; % min shaft shoulder diameter in mm
bear_thick=17; % bearing thickness in mm
spl= 300; %spool length in meters
g=9.81; %gravity in L/s^2
F_spl=50*g; % weight of spool and wire just before the winding is complete
F_pul=500;  % resultant pulley tension

Gear_ratio=0.6;
T1=30/pi; % Nm torque of the motor

T2=0.94*T1/Gear_ratio;

T_wire=T2/0.08; % 8 cm wire winding radius when spool is almost full


w= 76.5e-6;  % carbon steel(AISI 1050) unit weight in N/mm^3 (table A-5)
E=207*10^9; % carbon steel(AISI 1050) modulus of elasticity in Pa (table A-5)

if D_shaft>=41
    d_shaft=D_shaft;
    'Shoulder/shaft diameter set to default'
end

if Bearing_center_distance>bear_thick
    S=Bearing_center_distance;
else
    S=bear_thick; % minimum bearing center distance
    'Bearing distance set to default'
end


I_x=Inertia(d_shaft,S); % in m^4

W1= 0.25*w*pi*d_bearing^2;  % weight of shaft segment near the spool and ...
%pulley in N/mm

W2= 0.25*w*pi*d_shaft^2;  % weight of the shaft's middle segment N/mm


%all units in mm

H=bear_thick/2; %8.5 bearing thickness
I=spl/2; %150  distance to center of spool from end
J=spl+H; %308.5 distance to center of bearing 1 from end
K=spl+bear_thick; %317  distance to first shoulder
L=spl+S; %300+S distance to second shoulder
length= K+S+40;  % total length

%% Reaction Forces
[Ra_y,Rb_y,Ra_z,Rb_z]=Reaction_forces(W1,W2,S,T_wire,F_spl,F_pul);

%% Shear and Moment

V=@(x) -W1*ramp(x,0)-F_spl*step(x,I)+Ra_y*step(x,J)-(W2-W1)*ramp(x,K)-...
    (W1-W2)*ramp(x,L)+Rb_y*step(x,L+H)-F_pul*step(x,K+S+10);



M= @(x) -0.5*W1*parab(x,0)-F_spl*ramp(x,I)+Ra_y*ramp(x,J)- ...
    0.5*(W2-W1)*parab(x,K)-0.5*(W1-W2)*parab(x,L)+...
    Rb_y*ramp(x,L+H)-F_pul*ramp(x,K+S+10); % in N.mm



%theta_prime= @(x)(1/E)*M(x).*I_x(x).^-1;

I1=double(I_x(0));
I2=double(I_x(K+1));

syms C1 C0


% theta prime not including E

theta_prime=@(x) (1/I1)*(-0.5*W1*parab(x,0)-F_spl*ramp(x,I)+Ra_y*ramp(x,J))+...
    ((1/I2)-(1/I1))*(-0.5*W1*parab(x,K)-(K)*W1*ramp(x,K)-...
    0.5*W1*((K)^2)*step(x,K)-F_spl*ramp(x,K)-...
    (K-I)*F_spl*step(x,K)+Ra_y*ramp(x,K)+H*Ra_y*step(x,K))-...
    (W2-W1)*parab(x,K)/(2*I2); % in N.mm/m^4 just past first shoulder



theta_prime2=@(x) theta_prime(x)+((1/I1)-(1/I2))*(-0.5*W1*parab(x,L)-...
    L*W1*ramp(x,L)-0.5*W1*((L)^2)*step(x,L)-F_spl*ramp(x,L)-...
    (L-I)*F_spl*step(x,L)+Ra_y*ramp(x,L)+(S-H)*Ra_y*step(x,L)-...
    (W2-W1)*(0.5*parab(x,L)+(L-K)*ramp(x,L)+...
    0.5*(L-K)^2*step(x,L)))+(-0.5*(W1-W2)*parab(x,L)+...
    Rb_y*ramp(x,L+H)-F_pul*ramp(x,K+S+10))/I1;% in N.mm/m^4



Slope=@(x) (1/I1)*(-(1/6)*W1*Third(x,0)-0.5*F_spl*parab(x,I)+...
    0.5*Ra_y*parab(x,J))+((1/I2)-(1/I1))*(-(1/6)*W1*Third(x,K)-...
    0.5*(K)*W1*parab(x,K)-0.5*W1*((K)^2)*ramp(x,K)-...
    0.5*F_spl*parab(x,K)-(K-I)*F_spl*ramp(x,K)+0.5*Ra_y*parab(x,K)+...
    H*Ra_y*ramp(x,K))-(W2-W1)*Third(x,K)/(6*I2);% in N.mm^2/m^4

Slope2=@(x) Slope(x)+((1/I1)-(1/I2))*(-(1/6)*W1*Third(x,L)-...
    0.5*L*W1*parab(x,L)-0.5*W1*((L)^2)*ramp(x,L)-0.5*F_spl*parab(x,L)-...
    (L-I)*F_spl*ramp(x,L)+0.5*Ra_y*parab(x,L)+(S-H)*Ra_y*ramp(x,L)-...
    (W2-W1)*((1/6)*Third(x,L)+0.5*(L-K)*parab(x,L)+...
    0.5*(L-K)^2*ramp(x,L)))+(-(1/6)*(W1-W2)*Third(x,L)+...
    0.5*Rb_y*parab(x,L+H)-0.5*F_pul*parab(x,K+S+10))/I1;% in N.mm^2/m^4

Slope_var=@(x) (Slope2(x)/(E*10^6))+C1; %in rad

Defl=@(x) (1/I1)*(-(1/24)*W1*Fourth(x,0)-(1/6)*F_spl*Third(x,I)+...
    (1/6)*Ra_y*Third(x,J))+((1/I2)-(1/I1))*(-(1/24)*W1*Fourth(x,K)-...
    (1/6)*(K)*W1*Third(x,K)-0.5^2*W1*((K)^2)*parab(x,K)-...
    (1/6)*F_spl*Third(x,K)-0.5*(K-I)*F_spl*parab(x,K)+...
    (1/6)*Ra_y*Third(x,K)+0.5*H*Ra_y*parab(x,K))-...
    (W2-W1)*Fourth(x,K)/(24*I2); % in N.mm^3/m^4

Defl2=@(x) Defl(x)+((1/I1)-(1/I2))*(-(1/24)*W1*Fourth(x,L)-...
    (1/6)*L*W1*Third(x,L)-0.5^2*W1*((L)^2)*parab(x,L)-...
    (1/6)*F_spl*Third(x,L)-0.5*(L-I)*F_spl*parab(x,L)+...
    (1/6)*Ra_y*Third(x,L)+0.5*(S-H)*Ra_y*parab(x,L)-...
    (W2-W1)*((1/24)*Fourth(x,L)+(1/6)*(L-K)*Third(x,L)+...
    0.5^2*(L-K)^2*parab(x,L)))+(-(1/24)*(W1-W2)*Fourth(x,L)+...
    (1/6)*Rb_y*Third(x,L+H)-(1/6)*F_pul*Third(x,K+S+10))/I1;
% in N.mm^3/m^4

Defl_var=@(x) (Defl2(x)/(E*10^6))+C1*x+C0; %deflection in mm

def_b1=Defl_var(J);
def_b2=Defl_var(J+S);

Coeff=solve([def_b1==0,def_b2==0]);

C0_const=Coeff.C0;
C1_const=Coeff.C1;


Slope2=@(x) double((Slope2(x)/(E*10^6))+C1_const);%in rad

Defl2=@(x) double((Defl2(x)/(E*10^6))+C1_const*x+C0_const); %deflection in mm

end

