function [ I ] = Inertia( d_shaft, B_distance )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


bear_thick=17; %in mm

spl= 300; %spool length in mm

d_bearing=0.035; % in meters

d_shaft=d_shaft/1000; % change to meters

S=B_distance;


I1=(pi/64)*(d_bearing)^4;
I2=(pi/64)*(d_shaft)^4;

I= @(x) I1*step(x,0)+(I2-I1)*step(x,spl+bear_thick)+...
    (I1-I2)*step(x,spl+S);%assuming there is only one shoulder at each 
                          %bearing location & bearings are flat against 
                          %the spool in m^4
                          


end

