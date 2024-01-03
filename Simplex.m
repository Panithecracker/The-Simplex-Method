format long
%Insert any LP in SEF below:
%The program solves the LP, assuming the LP is not infeasible (2nd Phase
%Simplex)
A = zeros(5,7);
A(1,:) = [-1 1 1 0 0 0 0];
A(2,:) = [-0.5 1 0 1 0 0 0];
A(3,:) = [0.5 -1 0 0 1 0 0];
A(4,:) = [1 1 0 0 0 1 0];
A(5,:) = [1 0 0 0 0 0 1];
b = transpose([0 3 0 7 5.7]);
B = [3,4,5,6,7];
d = 0;
c = transpose([1 1 0 0 0 0 0]);
basis = [];
basic_sols = [];
[basis, basic_sols] = SimplexAlgo(c,d,A,b,B,basis,basic_sols);

%Visualization part:  (only applicable if the number of variables is <=3)
%Knowing the feasible basic solutions obtained in each iteration of the
%Simplex and knowing the polyhedral form of the LP, we can find the
%feasible solutions in the polyhedron related to those basic solutions at
%each iteration, by removing the slack from the basic solutions. This is
%more precisely a 1-1 mapping that allows to study the geometry of the
%simple algorithm, revealing many insights such as the fact that optimal
%LP's will always optimal solutions at the extreme points of the polyhedra and that Simplex  transformed basic solutions are always extreme points in the polyhedron.
%2D CASE:
P = zeros(5,2);
P(1,:) = [-1 1];
P(2,:) = [-0.5 1];
P(3,:) = [0.5 -1];
P(4,:) = [1 1];
P(5,:) = [1 0];
b = transpose([0 3 0 7 5.7]);
figure
xlabel("x1");
ylabel("x2");
title("Feasible region");
A=[1,-1;0.5,-1;-0.5,1;-1,-1];
b=transpose([0,-3,0,-7]);
lb = transpose([-inf -inf]);
rb = transpose([5.7 inf]);
plotregion(A,b,lb,rb,[0.7,0.2,0.3]);
axis equal
grid on;
hold on;
h = scatter(basic_sols(1,1),basic_sols(2,1),"red","*");
hold on;
for j=2:size(basic_sols,2)
    pause(0.5);
    h.XData = basic_sols(1,j);
    h.YData = basic_sols(2,j);
    hold on;
end

function [mem_b,mem_s] = SimplexAlgo(c,d,A,b,B, mem_b, mem_s)
%% input details:
%LP in SEF given by:
%1- objective function : made up by a vector c and a number d
%2- matrix A and vector b
%3- feasible basis (vector of column indices in increasing order )
%4- mem_b and mem_s are matrices whose columns are the feasible basis and respective solutions used by simplex
%note: mem_b has a varying amount of columns but the dimension of the rows
%are the number of rows. For mem_s the same happens but the dimension of rows
%equals the number of columns.
%until terminating
%% output details:
% type : a string with values "unbounded" or "optimal"
% certificate : a vector with respective certificate data
c_b = zeros(size(B,2),1);
c_n = zeros(size(A,2) - size(B,2),1);
%% Step 1: Convert the LP into Canonical Form for basis B
% + efficiency note : check if it is in CF for B already...
A_b = zeros(size(A,1),size(A,1));
for i=1:size(A,1) %get the submatrix for the basis B
    A_b(:,i) = A(:,B(i)); 
end
A_binv = inv(A_b);
%get the associated objective function for the new reg
[c_b, c_n] = sub_divide(c,B); 
y = (transpose(inv(A_b)))*c_b;
d = d + transpose(y)*b;
c = transpose(transpose(c) - (transpose(y)*A));
%multiply both sides of the system Ax = b by its inverse for constraint reg
A = A_binv*A;
b = A_binv*b;
% After this, the LP is in canonical form for feasible B.
fprintf("Basis : \n");
disp(B);
mem_b(:,size(mem_b,2)+1) = transpose(B); %update the current basic solution
fprintf("Basic solution : \n");
bsol = zeros(size(A,2),1); %basic solution of current basis B
% 
 for i=1:size(B,2) %we get the basic solution for now
    bsol(B(1,i),1) = b(i,1);
 end
mem_s(:,size(mem_s,2)+1) = bsol;
disp(bsol)
fprintf("Value : %d\n", transpose(c)*bsol+d);
% 
 %% Step 2: Optimality test
 %subdivide c into c_b and c_n for later usage:
[c_b , c_n] = sub_divide(c,B);
positive = 0; %number of positive entries for current c
for i = 1:size(c_n,1) %check if c_n <= 0 (c_b is zero by definition of Cform for B)
    if(c_n(i,1) > 0)
        positive = positive +1;
    end
end
if (positive == 0) %then c <= 0 so we have found an optimal solution
    disp("The LP is optimal");
else
    % the current basic solution is not optimal at this point
    %-----
    %% Step 3: Use blands rule and test for unboundedness before anything
    for i=1:size(c,1)
        if (~ismember(i,transpose(B)))&&(c(i) > 0)
            k = i; %then the ith variable will enter the basis... 
            break
        end
    end
    A_k = A(:,k); %get the kth column of A 
    if( sum(A_k <= 0,1) == size(A_k,1)) %if A_k <= 0 then stop! (the LP is unbounded)
        disp("The LP is unbounded with certificate : ");
    else
        % up to this point, the optimality and unboundedness tests have failed.
        %Then, we can make another iteration of Simplex , finding a new
        %basis whose basic solution will have a value greater or equal than
        %that of the current basic solution. The equal case is called a
        %degenerate case. Also, when determining the entering and exitting variables we
        %use Bland's Rule which is proven to avoid cycling.
        min_ratio = inf;
        for i=1:size(A_k,1)
            if(A_k(i,1) > 0)
                ratio = (b(i,1))/(A_k(i,1));
                if (ratio < min_ratio)
                    r = i;
                    min_ratio = ratio;
                end
            end
        end
        fprintf("k = %d enters basis, l = %d leaves basis\n",k,B(1,r));
        %update basis for next iteration:
        B(1,r) = k; 
        B = sort(B);
        %Next iteration
        [mem_b, mem_s]=SimplexAlgo(c,d,A,b,B, mem_b, mem_s);
     end
 end
end


%% useful complementary functions

function [v_b, v_n ] = sub_divide(v, B) %given a column and row of indices B, returns the subdivision into those in B and not in B
   v_b = zeros(size(B,2),1);
   v_n = zeros(size(v,1) - size(B,2),1);
   temp = v;
   for i=1:size(B,2)
       v_b(i,1) = v(B(1,i)); 
       temp(B(1,i)) = inf; %mark those entries in b from the input vector for later
   end
   j = 1;
   pos = 1;
   while (pos <= size(v_n,1))
       if(~isinf(temp(j,1)))
           v_n(pos,1) = temp(j,1);
           pos = pos +1;
       end
       j = j +1;
   end
end



