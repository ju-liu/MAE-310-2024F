clear all; clc; % clean the memory and screen

% Define the external source or force and boundary data
f = @(x) x; % f(x) = x
g = 1.0;    % u    = g  at x = 1
h = 0.5;    % -u,x = h  at x = 0

% Setup the mesh
pp   = 1;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_el = 4;              % number of elements
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations
n_int = pp;

hh = 1.0 / (n_np - 1); % space between two adjacent nodes
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

IEN = zeros(n_el, n_en);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(ee, aa) = (ee - 1) * pp + aa;
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% allocate the stiffness matrix
K = zeros(n_eq, n_eq);

for ee = 1 : n_el
  
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
  
  %x_ele = zeros(n_en, 1);
  %for aa = 1 : n_en
  %  x_ele(aa) = x_coor(IEN(ee,aa));
  %end
  x_ele = x_coor(IEN(ee,:));

  for qua = 1 : n_int
    
    dx_dxi = 0.0;
    for aa = 1 : n_en
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
    end
    dxi_dx = 1.0 / dx_dxi;

    for aa = 1 : n_en
      for bb = 1 : n_en
        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
      end
    end
  end
  
  % check the ID(IEN(ee, aa)) and ID(IEN(ee,bb)), if they are positive
  % put the element stiffness matrix into K
  for aa = 1 : n_en
    P = ID(IEN(ee,aa));
    if(P > 0)
      for bb = 1 : n_en
        Q = ID(IEN(ee,bb));
        if(Q > 0)
          K(P, Q) = K(P, Q) + k_ele(aa, bb);
        end
      end
    end
  end

end















% EOF