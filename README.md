# The-Simplex-Method
Implementation of Simplex Method, used to find the solution of any type of optimization problem, where the objective function is affine and the constraints are linear. 
# Visual Example:
In order to illustrate the geoemetrical meaning of the algorithm while showing the general functionality of this project, I applied the Simplex Algorithm to an LP whose feasible region is a regular polygon with n sides. As the gif below shows, each iteration of the simplex method computes a new feasible solution which is not coincidentially, an extreme point of the polyhedron (feasible region for the example). This is in fact always the case: Basic feasible solutions are extreme points of the polyhedra and any optimal LP attains its optimal value at some extreme point(s). For this case, the objective function is $f(x,y) = x+y$ and the program succesfully returns the optimal solution accordingly for the given feasible region (the labeled points are the solutions that Simplex produced; the last one was $(240.45,129.39)$ which is optimal).
Note: The code for representing this feasible region in matrix form is also uploaded in the repository under the name "VisualExample.m"

![image](https://github.com/Panithecracker/The-Simplex-Method/assets/97905110/a9361e74-f111-4453-bd12-e8b10f206eac)


