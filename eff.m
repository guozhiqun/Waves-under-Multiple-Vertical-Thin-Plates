function [y] = eff(x,h,jp0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y=x*tan(x*h)+(jp0)^2/9.80665;
end


