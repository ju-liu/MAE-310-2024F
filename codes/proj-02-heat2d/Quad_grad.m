function [val_xi, val_eta] = Quad_grad(aa, xi, eta)

if aa == 1
    val_xi  = -0.25 * (1-eta);
    val_eta = -0.25 * (1-xi);
elseif aa == 2
    val_xi  =  0.25 * (1-eta);
    val_eta = -0.25 * (1+xi);
elseif aa == 3
    val_xi  = 0.25 * (1+eta);
    val_eta = 0.25 * (1+xi);
elseif aa == 4
    val_xi  = -0.25 * (1+eta);
    val_eta =  0.25 * (1-xi);
else
    error('Error: value of a should be 1,2,3, or 4.');
end

% EOF