function [ theta_prime t_prime2 ] = Top_Shaft_test( D_shaft, Bearing_center_distance )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

d_bearing=35; % bearing bore diameter in mm
d_shaft=45; % min shaft shoulder diameter in mm
bear_thick=17; % bearing thickness in mm
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

if D_shaft>=45
    d_shaft=D_shaft;
end

if Bearing_center_distance>bear_thick
    S=Bearing_center_distance;
else
    S=bear_thick; % minimum bearing center distance
end


I_x=Inertia(d_shaft,S);

W1= w*pi*d_bearing^2;  % weight of shaft segment near the spool and ...
                            %pulley in N/mm
                            
W2= w*pi*d_shaft^2;  % weight of the shaft's middle segment 



[Ra_y,Rb_y,Ra_z,Rb_z]=Reaction_forces(W1,W2,S,T_wire,F_spl,F_pul);


H=bear_thick/2; %8.5
I=spl/2; %150
J=spl+H; %308.5
K=spl+bear_thick; %317
L=spl+S; %300+S
length= K+S+40;

V=@(x) -W1*ramp(x,0)-F_spl*step(x,I)+Ra_y*step(x,J)-(W2-W1)*ramp(x,K)-...
       (W1-W2)*ramp(x,L)+Rb_y*step(x,L+H)-F_pul*step(x,K+S+10);


    
M= @(x) -0.5*W1*parab(x,0)-F_spl*ramp(x,I)+Ra_y*ramp(x,J)- ...
        0.5*(W2-W1)*parab(x,K)-0.5*(W1-W2)*parab(x,L)+...
        Rb_y*ramp(x,L+H)-F_pul*ramp(x,K+S+10);
   

    
theta_prime= @(x)(1/E)*M(x).*I_x(x).^-1;

I1=double(I_x(0));
I2=double(I_x(318));



t_prime=@(x) (1/I1)*(-0.5*W1*parab(x,0)-F_spl*ramp(x,I)+Ra_y*ramp(x,J))+...
        ((1/I2)-(1/I1))*(-0.5*W1*parab(x,K)-(K)*W1*ramp(x,K)-...
        0.5*W1*((K)^2)*step(x,K)-F_spl*ramp(x,K)-...
        (K-I)*F_spl*step(x,K)+Ra_y*ramp(x,K)+H*Ra_y*step(x,K))-...
         (W2-W1)*parab(x,K)/(2*I2);
     


t_prime1=@(x) t_prime(x)+((1/I1)-(1/I2))*(-0.5*W1*parab(x,L)-...
           L*W1*ramp(x,L)-0.5*W1*((L)^2)*step(x,L)-F_spl*ramp(x,L)-...
         (L-I)*F_spl*step(x,L)+Ra_y*ramp(x,L)+(S-H)*Ra_y*step(x,L)-...
         (W2-W1)*(0.5*parab(x,L)+(L-K)*ramp(x,L)+...
         0.5*(L-K)^2*step(x,L)))+(-0.5*(W1-W2)*parab(x,L)+...
         Rb_y*ramp(x,L+H)-F_pul*ramp(x,K+S+10))/I1;
        

t_prime2=@(x) t_prime1(x)/E;
    
end

