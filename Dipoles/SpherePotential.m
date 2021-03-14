function [F,X,Y] = SpherePotential(XYZ,Q,R,r0,a,b,Dx,Dy,Nxy)
    n = length(Q);
    e1 = a / ((a(1,1)^2+a(2,1)^2+a(3,1)^2)^0.5);
    c2 = zeros(3,1);
    c2(1, 1) = b(2,1) * a(3,1) - a(2,1) * b(3,1);
    c2(2, 1) = b(3,1) * a(1,1) - a(3,1) * b(1,1);
    c2(3, 1) = b(1,1) * a(2,1) - a(1,1) * b(2,1);
    e2 = c2 / ((c2(1,1)^2+c2(2,1)^2+c2(3,1)^2)^0.5);
    X = zeros(Nxy(1),1);
    Y = zeros(Nxy(2),1);
    F = zeros(Nxy(1),Nxy(2));
    for ii = 1:Nxy(1)
        X(ii, 1) = Dx(1) + ((Dx(2) - Dx(1)) / Nxy(1))*ii;
    end
    for jj = 1:Nxy(2)
        Y(jj, 1) = Dy(1) + ((Dy(2) - Dy(1)) / Nxy(2))*jj;
    end
    for ii = 1:Nxy(1)
        for jj = 1:Nxy(2)
            for ll = 1:n
                if (((r0(1,1) + X(ii, 1) * e1(1,1) + Y(jj, 1) * e2(1,1) - XYZ(ll,1))^2 + (r0(2,1) + X(ii, 1) * e1(2,1) + Y(jj, 1) * e2(2,1) - XYZ(ll,2))^2 + (r0(3,1) + X(ii, 1) * e1(3,1) + Y(jj, 1) * e2(3,1) - XYZ(ll,3))^2)^(0.5) < 1)
                    F(ii, jj) =  F(ii, jj) + Q(ll,1) / R(ll);
                else
                    F(ii, jj) =  F(ii, jj) + Q(ll,1) / ((r0(1,1) + X(ii, 1) * e1(1,1) + Y(jj, 1) * e2(1,1) - XYZ(ll,1))^2 + (r0(2,1) + X(ii, 1) * e1(2,1) + Y(jj, 1) * e2(2,1) - XYZ(ll,2))^2 + (r0(3,1) + X(ii, 1) * e1(3,1) + Y(jj, 1) * e2(3,1) - XYZ(ll,3))^2)^(0.5);
                end
            end
        end
    end    
end
