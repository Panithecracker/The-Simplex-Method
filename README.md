# The-Simplex-Method
Implemmentation of Simplex Method, used to find the solution of any type of optimization problem, where the objective function is affine and the constraints are linear. 
# Visual Example:
In order to illustrate the geoemetrical meaning of the algorithm while showing the general functionality of this project, I applied the Simplex Algorithm to an LP whose feasible region is a regular polygon with n sides. As the gif below shows, each iteration of the simplex method computes a new feasible solution which is not coincidentially, an extreme point of the polyhedron (feasible region for the example). This is in fact always the case: Basic feasible solutions are extreme points of the polyhedra and any optimal LP attains its optimal value at some extreme point(s). For this case, the objective function is f(x,y) = x+y and the program succesfully returns the optimal solution accordingly (it is the last extreme point colored in blue in the animation).

