format long
%% Example : LP with n-regular polygon in 2D as a feasible region :
%system of linear constraints in polyhedral form.
% For the example, the feasible region is a n-gon centered at (200,100) with
% pivot at % (250,100).
n = 10; 
V = CreateNgon(n,200,100,250,100); %vertices of the n-gon
[P,b] = FindConstraints(n,V); %now, the feasible region is Px <= b
figure
scatter(V(1,:)',V(2,:)',"blue","filled");
xlabel("x1");
ylabel("x2");
title("Feasible region plot");
hold on;
% 
% %testing that the equations are correct for the feasible region
% x = linspace(150,250,1000);
% y = zeros(1,length(x));
% for l=1:n
%     for i=1:1000
%         y(i) = (-1)*((P(l,1))/(P(l,2)))*x(i) + (b(l))/(P(l,2));
%     end
%     plot(x,y, "blue");
%     y = [];
%     hold on;
% end
% axis equal



plotregion(-P,-b,[0 0], [inf inf], [0.8 1 0.9]); %times -1 to express it as -Px >=-b for the format of plotregion command
grid on;
hold on;

% Convert the LP into SEF and use Simplex function to find the solution
A = [P eye(n)];
B = [3:n+2];
d = 0;
c = zeros(n+2,1);
c(1,1) = 1; %objective function f(x,y,s1,...,sn) = x+y;
c(2,1) = 1;
basis = [];
basic_sols = [];
%finally make sure the b vector is nonnegative 
for r=1:size(A,1)
    if(b(r,1) < 0)
        A(r,:) = (-1)*A(r,:);
        b(r,1) = (-1)*b(r,1);
    end
end
disp("In the beginning, the matrix A is : ");
A
[basis, basic_sols] = SimplexAlgo(c,d,A,b,B,basis,basic_sols);

%Animate it
h = scatter(basic_sols(1,1),basic_sols(2,1),[],[0.1 0 0.8], "filled");
hold on;
caption = sprintf('x1=%.2f, x2=%.2f',basic_sols(1,1) , basic_sols(2,1));
text(basic_sols(1,1)+0.05, basic_sols(2,1)+0.05, caption);
for j=2:size(basic_sols,2)
    pause(1);
    h.XData = basic_sols(1,j);
    h.YData = basic_sols(2,j);
    caption = sprintf('x1=%.2f, x2=%.2f',basic_sols(1,j) , basic_sols(2,j));
    text(basic_sols(1,j)+0.05, basic_sols(2,j)+0.05, caption);
    hold on;
end

