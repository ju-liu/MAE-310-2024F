clear all; clc;

kappa = 1.0 % conductivity

exact = @(x,y) x*(1-x)*y*(1-y); % exact solution
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) 2*kappa*x*(1-x)+2*kappa*y*(1-y); % source
