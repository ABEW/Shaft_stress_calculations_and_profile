function [ par ] = parab(x,a)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if(x>a)
    par=(x-a)^2;
else
    par=0;
end

