function [extreme_points] = CreateNgon(n,c1,c2, p1,p2) 
%Given the number of sides, the center and a pivot , this function returns
%a matrix whose columns are the extreme points of the respective ngon.
%Step1 : Find the different vertices of a regular ngon with n sides of
%length l centered at the origin
theta = (2*pi)/(n);
%use the rotation matrix by theta radians iteratively as follows:
R(1,:) = [cos(theta) -sin(theta)];
R(2,:) = [sin(theta) cos(theta)];
extreme_points = zeros(2,n);
extreme_point = transpose([p1-c1, p2-c2]); %first extreme point is the given pivot (usually center + [l 0]')
for i=1:n
    extreme_points(1,i) = extreme_point(1,1);
    extreme_points(2,i) = extreme_point(2,1);
    extreme_point = R*extreme_point;
end
%Now, retranslate it to the center given to obtain the requested polygon
for j=1:n
    extreme_points(1,j) = extreme_points(1,j) + c1;
    extreme_points(2,j) = extreme_points(2,j) +c2;
end
end
