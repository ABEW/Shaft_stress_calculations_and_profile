function [ S ] = step2( x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

del=1e-16;

S=heaviside(x+del);
end

