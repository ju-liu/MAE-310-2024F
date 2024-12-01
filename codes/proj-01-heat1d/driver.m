clear all; clc; clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

% Setup the mesh
pp   = 2;              % polynomial degree
n_en = pp + 1;         % number of element or local nodes
n_el = 2;              % number of elements
n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations
n_int = 10;

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
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
F = zeros(n_eq, 1);

% Assembly of the stiffness matrix and load vector
for ee = 1 : n_el
    k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
    f_ele = zeros(n_en, 1);    % allocate a zero element load vector

    x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

    % quadrature loop
    for qua = 1 : n_int
        dx_dxi = 0.0;
        x_l = 0.0;
        for aa = 1 : n_en
            x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
        end
        dxi_dx = 1.0 / dx_dxi;

        for aa = 1 : n_en
            f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
            for bb = 1 : n_en
                k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
            end
        end
    end

    % Assembly of the matrix and vector based on the ID or LM data
    for aa = 1 : n_en
        P = ID(IEN(ee,aa));
        if(P > 0)
            F(P) = F(P) + f_ele(aa);
            for bb = 1 : n_en
                Q = ID(IEN(ee,bb));
                if(Q > 0)
                    K(P, Q) = K(P, Q) + k_ele(aa, bb);
                else
                    F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
                end
            end
        end
    end
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

% Solve Kd = F equation
d_temp = K \ F;

disp = [d_temp; g];

% Postprocessing: visualization
%plot(x_coor, disp, '--r','LineWidth',3);

%x_sam = 0 : 0.01 : 1;
%y_sam = x_sam.^5;
%hold on;
%plot(x_sam, y_sam, '-k', 'LineWidth', 3);

n_sam = 20;
xi_sam = -1 : (2/n_sam) : 1;

x_sam = zeros(n_el * n_sam + 1, 1);
y_sam = x_sam; % store the exact solution value at sampling points
u_sam = x_sam; % store the numerical solution value at sampling pts

for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, :) );
    u_ele = disp( IEN(ee, :) );

    if ee == n_el
        n_sam_end = n_sam+1;
    else
        n_sam_end = n_sam;
    end

    for ll = 1 : n_sam_end
        x_l = 0.0;
        u_l = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
            u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
        end

        x_sam( (ee-1)*n_sam + ll ) = x_l;
        u_sam( (ee-1)*n_sam + ll ) = u_l;
        y_sam( (ee-1)*n_sam + ll ) = x_l^5;
    end
end


plot(x_sam, u_sam, '-r','LineWidth',3);
hold on;
plot(x_sam, y_sam, '-k','LineWidth',3);


% el2 = sqrt(integral(fun,0,1)) / sqrt(integral(fun2,0,1))
% x = sym ('x');
% function f = @(x) -20*x.^3
% el2 = (integral(x,0,1)) / (integral(x,0,1))
% x = sym("x");
% y_sam 是 理论解 黑色
% u_sam 是 近似解 红色

% % u_sam =
% % y_sam =
% 
% x = u_sam - y_sam;
% f = @(x) (x)^2;
% el2 = (integral(f,0,1)) / (integral(f,0,1))


% x_error = u_sam - y_sam;
% f1 = @(x_error) (x_error).^2;
% f2 = @(y_sam) (y_sam).^2;
% el2 = sqrt(integral(f1,0,1)) / sqrt(integral(f2,0,1))

% x_error = u_sam - y_sam;
% f1 = @(x_error) (x_error).^2;
% f2 = @(y_sam) (y_sam).^2;
% el2 = sqrt(integral(f1(x_error),0,1)) / sqrt(integral(f2(y_sam),0,1))

% x_error = u_sam - y_sam;
% f1 = @(x) (x).^2;
% el2 = sqrt(integral(@(x) f1(x_error),0,1)) / sqrt(integral(@(x) f1(y_sam),0,1))

% x_error = u_sam - y_sam;
% f1 = @(x) (x).^2;
% el2 = sqrt(integral(f1(x_error),0,1)) / sqrt(integral(f1(y_sam),0,1))


% % el2 = sqrt(integral(fun,0,1)) / sqrt(integral(fun2,0,1))

% fun = @(x,c) 1./(x.^3-2*x-c);
% q = integral(@(x) fun(x,5),0,2)
% eh1

sum = 0;
A = [1 2 3];
for counting = 1 : 3
    f = @(x) A(counting).^2;
    sum = integral(f,0,1) + sum;
end
sum


% EOF