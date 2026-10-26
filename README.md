# Design-and-Analysis-of-the-Capacitance-of-an-Arbitrarily-Shaped-Capacitor

This project presents a computational tool for the design and analysis of the capacitance of an arbitrarily shaped capacitor. Many numerical methods can be implemented to solve this
problem, and each has their respective strengths and limitations for different cases. With arbitrarily shaped capacitors and parasitic capacitance effects in modern analog electronics, it is useful to understand under which circumstances would demand implementing a certain numerical method. We built and evaluated numerical solvers of the Laplace Equation using two different numerical methods: the Finite Difference Method (FDM) and the Boundary Element Method (BEM). The objective of this project is to determine the merits and limitations of each numerical method for different geometries and circumstances, and to provide insight into how these models could be used in analog circuit design and layout.

## Finite Difference Method (FDM)
The Finite Difference Method (FDM) is a numerical technique for approximating solutions to differential equations that involves discretizing the solution domain into a grid, and
setting up a system of linear equations relating the derivatives at different grid points to one another using Taylor series approximations [1]. 

Our implementation of the FDM solver is found in ```LaplaceSolver3DSurfacesFDM.m```. To run the solver on your own device, download the respective ```.mlx``` file from the ```Final_Submission``` folder. For ease of testing, the ```FDMTestFiles``` folder has example test files for different surfaces. To test different results seen in our report, download the respective test file in MATLAB and set the different conditions accordingly. More details can be found in the comments at the top of each test file.

## Boundary Element Method (BEM)
The Boundary Element Method (BEM) is a numerical method for finding the solution to a PDE, often used as an alternative to the Finite Element Method, or FEM, when the FEM solution would require an unreasonably large number of mesh points. BEM reduces the number of mesh points required by evaluating the PDE only on the boundaries of its domain, while FEM evaluates the PDE at mesh points throughout the entire volume [2].

Our implementation of the BEM solver is found in ```MAT292BEMsolver.m```. To run the solver on your own device, download the respective ```.mlx``` file from the ```Final_Submission``` folder. To test the different results seen in our report, please set z_bot and z_top fields, located on lines 47 and 48, equal to the case you would like to test.

For example, to test a parallel plate capacitor, you should have:

z_bot = @(x, y) z_bot_flat (x, y);

z_top = @(x, y) z_top_flat (x, y); 

You can change the length, width, and height using the Lx, Ly, and Lz fields. 



## Citations
[1] S. Chowdhury, P. K. Das, and S. B. Faruque, Numerical Solutions of Boundary Value Problems with Finite Difference Method. San Rafael, California: Morgan & Claypool Publishers; Institute of Physics, 2018.

[2] T. LaForce, “Boundary element method: Course notes,” Stanford University, Stanford, CA, Tech. Rep., Jun. 2006, pE281 Boundary Element Method course notes
