function [Q,D] = ElectroStaticDipoles(XYZ,R,F)
n=length(F);
A = zeros(4*n,4*n);
for ii=1:n
    for jj=1:n
        if (ii==jj)
            A(ii,jj)=1/R(ii);
        else
            A(ii,jj)=((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-0.5);
        end
        
    end
end

for ii=1:n
    for jj=1:n
        if (ii~=jj)
            A(ii,n+(jj-1)*3+1)=(XYZ(ii,1)-XYZ(jj,1))*((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-1.5);
            A(ii,n+(jj-1)*3+2)=(XYZ(ii,2)-XYZ(jj,2))*((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-1.5);
            A(ii,n+(jj-1)*3+3)=(XYZ(ii,3)-XYZ(jj,3))*((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-1.5);
        end
    end
end


for ii=1:n
    for jj=1:n
        if (ii~=jj)
            A(n+(jj-1)*3+1,jj)=R(ii)^3*(XYZ(ii,1)-XYZ(jj,1))*((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-1.5);
            A(n+(jj-1)*3+2,jj)=R(ii)^3*(XYZ(ii,2)-XYZ(jj,2))*((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-1.5);
            A(n+(jj-1)*3+3,jj)=R(ii)^3*(XYZ(ii,3)-XYZ(jj,3))*((XYZ(ii,1)-XYZ(jj,1))^2+(XYZ(ii,2)-XYZ(jj,2))^2+(XYZ(ii,3)-XYZ(jj,3))^2)^(-1.5);
        end
    end
end
for ii=n+1:4*n
    A(ii,ii)=-1;
end
    
f=zeros(4*n,1);
for k=1:n
    f(k)=F(k);
end
X=inv(A)*f;
Q=zeros(n,1);
D=zeros(3*n,1);
for k=1:n
    Q(k)=X(k);
end
for k=1:3*n
    D(k)=X(k+n);
end
end