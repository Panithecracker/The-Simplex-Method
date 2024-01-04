function [A,b] = FindConstraints(n,P)
%Given the n extreme points of a regular ngon in 2D, this function returns the
%n*2 matrix A and n*1 column vector b to represent the interior of the ngon
%in polyhedral form ie: Ax<= b 
%
A = zeros(n,2);
b = zeros(n,1);
for k=1:n
    if (k <=n/2)
        R(k) = (P(2,k)-P(2,k+1))/(P(1,k)-P(1,k+1)); %slope of kth side
        A(k,1) = R(k);
        A(k,2) = -1;
        b(k,1) = R(k)*P(1,k+1) - P(2,k+1);
    else
        if k ~= n
            R(k) = (P(2,k)-P(2,k+1))/(P(1,k)-P(1,k+1)); %slope of kth side
            A(k,1) = (-1)*R(k);
            A(k,2) = 1;
            b(k,1) = P(2,k+1) -R(k)*P(1,k+1);
        else
            R(k) = (P(2,1)-P(2,k))/(P(1,1)-P(1,k));
            A(k,1) = (-1)*R(k);
            A(k,2) = 1;
            b(k,1) = P(2,1) -R(k)*P(1,1);
        end
    end
    
end
A = -A;
b = -b;
end
