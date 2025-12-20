This project presents a computational tool for the design and analysis of the capacitance of an arbitrarily shaped capacitor. 
We built and tested numerical solvers using two different numerical methods: the Finite Difference Method (FDM) and the Boundary Element Method (BEM). 
Both of these solvers are availiable here! You will find them in the folder titled Final Submission. To run them, please download and open the appropriate file in MATLAB. Feel free to play around with the visualizations!

For the BEM solver, there is only one file you need to download. To test the different results seen in our report, please set z_bot and z_top fields, located on lines 47 and 48, equal to the case you would like to test.
For example, to test a parallel plate capacitor, you should have:

z_bot = @(x, y) z_bot_flat (x, y);

z_top = @(x, y) z_top_flat (x, y); 

You can change the length, width, and height using the Lx, Ly, and Lz fields. 

For the FDM solver, there is one main file you need to download. For ease of testing, FDMTestFiles folder has example test files for different surfaces. To test different results seen in our report, refer to the respective test file and set the different conditions accordingly. More details can be found in the comments at the top of each test file.
