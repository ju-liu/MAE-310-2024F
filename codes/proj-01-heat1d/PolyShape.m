function val = PolyShape(degree, a, xi, der)
switch degree
    % linear basis function
    case 1
        if a == 1
            if der == 0
                val = 0.5 * (1-xi);
            elseif der == 1
                val = -0.5;
            end
        elseif a == 2
            if der == 0
                val = 0.5 * (1+xi);
            elseif der == 1
                val = 0.5;
            end
        end

    % quadratic basis function
    case 2
        if a == 1
            if der == 0
                val = 0.5 * xi * (xi-1);
            elseif der == 1
                val = xi-0.5;
            end
        elseif a == 2
            if der == 0
                val = 1 -xi^2;
            elseif der == 1
                val = -2 * xi;
            end
        elseif a == 3
            if der == 0
                val = 0.5 * xi * (xi+1);
            elseif der == 1
                val = xi + 0.5;
            end
        end

    % cubic basis function
    case 3
        if a == 1
            if der == 0
                val = -9*(xi-(1/3)) * (xi+(1/3)) * (xi-1)/16;
            elseif der == 1
                val = -9*(2*xi*(xi-1)+xi^2-(1/9))/16;
            end
        elseif a == 2
            if der == 0
                val = 27 *(xi^2-1)*(xi-(1/3))/16;
            elseif der == 1
                val = 27 * (2*xi*(xi-(1/3))+xi^2-1)/16;
            end
        elseif a == 3
            if der == 0
                val = -27 * (xi^2-1)*(xi+(1/3))/16;
            elseif der == 1
                val = -27 * (2*xi*(xi+(1/3))+xi^2-1)/16;
            end
        elseif a == 4
            if der == 0
                val = 9*(xi+1)*(xi^2-(1/9))/16;
            elseif der == 1
                val = 9*(xi^2-(1/9)+2*xi*(xi+1))/16;
            end
        end

    otherwise
        error('Error: degree has to be 1, 2, or 3.');
end
