function [Q] =  ElectroStaticBalls(XYZ,R,F)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n=length(F);
A = ones(n,n);
for ii=1:n
    for jj=1:n
        if (ii==jj)
            A(ii,jj)=1/R(ii);
        else
            A(ii,jj)=((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-0.5);
        end
    end
end
Q=inv(A)*F;
end
