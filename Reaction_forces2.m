function [Ra_y,Rb_y,Ra_z,Rb_z]=Reaction_forces2(W1,W2,W3,W4,db1,db2,dfpul,...
                                    dfspl,K,L,N,O,T_wire,F_spl,F_pul)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Rb_y = (-0.5*W1*(db1)^2-F_spl*(db1-dfspl)+0.5*W1*(K-db1)^2+...
    W2*(L-K)*(0.5*(L+K)-db1)+W3*(N-L)*(0.5*(N+L)-db1)+F_pul*(dfpul-db1)+...
    W4*((O-N)*(0.5*(O+N)-db1)))*(db2-db1)^-1;
Ra_y = W1*K+W2*(L-K)+W3*(N-L)+W4*(O-N)+F_pul+F_spl-Rb_y;


Rb_z=-T_wire*db1/(db2-db1);

Ra_z=T_wire-Rb_z;

end

