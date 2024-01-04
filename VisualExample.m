format long
%% Example : LP with n-regular polygon in 2D as a feasible region :
%system of linear constraints in polyhedral form.
% For the example, the feasible region is a n-gon centered at (40,30) with
% pivot at % (50,30).
n = 10; 
V = CreateNgon(n,40,30,50,30); %vertices of the n-gon
[P,b] = FindConstraints(n,V); %now, the feasible region is Px <= b
A = [P eye(n)];
B = [3:n+2];
d = 0;
c = zeros(n+2,1);
c(1,1) = 1;
c(2,1) = 1; 
basis = [];
basic_sols = [];
[basis, basic_sols] = SimplexAlgo(c,d,A,b,B,basis,basic_sols);
figure
scatter(V(1,:)',V(2,:)',"white","filled");
xlabel("x1");
ylabel("x2");
title("Feasible region");
plotregion(-P,-b,[-inf -inf ], [inf inf], [0.8 1 0.9]); %times -1 to express it as -Px >=-b for the format of plotregion command
grid on;
hold on;
h = scatter(basic_sols(1,1),basic_sols(2,1),[],[0.1 0 0.8], "filled");
hold on;
for j=2:size(basic_sols,2)
    pause(0.5);
    h.XData = basic_sols(1,j);
    h.YData = basic_sols(2,j);
    hold on;
end

