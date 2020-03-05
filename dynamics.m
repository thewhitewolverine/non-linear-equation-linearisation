function [ xdot ] = dynamics(t, z)
%DYNAMICS Summary of this function goes here
%   Detailed explanation goes here
% Satya 16D170026   Exam 24
global U_k

h(1) = z(1);
h(2) = z(2);
h(3) = z(3);

U(1) = U_k(1);
U(2) = U_k(2);

hdot(1) = -U(1)*h(1) + h(1)*((0.48*h(2)*(1-h(3)/50))/(1.2+h(2)+(h(2)^2)/22));
hdot(2) = U(1)*(U(2)-h(2)) - 2.5*((0.48*h(2)*(1-h(3)/50))/(1.2+h(2)+(h(2)^2)/22));
hdot(3) = -U(1)*h(3) + h(1)*(0.2 + 2.2*((0.48*h(2)*(1-h(3)/50))/(1.2+h(2)+(h(2)^2)/22)));

xdot = hdot';