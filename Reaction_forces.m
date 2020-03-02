function [Ra_y Rb_y Ra_z Rb_z ]=Reaction_forces(w1,w2,S,T_wire,F_spl,F_pul )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Rb_y= -(1/S)*(w1*(0.5*(308.5^2)-57*(S+20)-0.5*8.5^2)-w2*S*(S-17)/2+...
    F_spl*158.5-F_pul*(S+18.5));

Ra_y=374*w1+w2*(S-17)+F_spl+F_pul-Rb_y;


Rb_z= -T_wire*(158.5)/S;

Ra_z= T_wire-Rb_z;

end

