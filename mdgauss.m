function p = mdgauss(x,mu,C)

[m,n] = size(x);
iC = inv(C+eps);
dC = det(C+eps);
k = 1/sqrt((2*pi)^m*dC);

p=zeros(n,1);

for i=1:n
    p(i)=k*exp(-1/2*(x(:,i)-mu)'*iC*(x(:,i)-mu));
end

end


