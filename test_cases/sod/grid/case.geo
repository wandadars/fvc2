// This GMSH file is drawn in units of meters. It is a 1D grid for a SOD shock tube case.

tube_length = 1; // Length of the shock tube
tube_height = 0.1; // Height of the shock tube (irrelevant, but needed for 2D geometry)
num_cells = 10000; // Number of cells along the length of the tube

// Define points starting in lower left corner and going clockwise around the geometry
Point(1) = {0.0, 0.0, 0.0}; // Bottom left corner
Point(2) = {0.0, tube_height, 0.0}; // Top left corner
Point(3) = {tube_length, tube_height, 0.0}; // Top right corner
Point(4) = {tube_length, 0.0, 0.0}; // Bottom right corner


// Define lines that make up the domain
Line(1) = {1, 2}; // Inlet
Line(2) = {2, 3}; // Top wall
Line(3) = {3, 4}; // Outlet
Line(4) = {4, 1}; // Bottom wall


Curve Loop(1) = {1, 2, 3, 4};
Surface(1) = {1};
Transfinite Line {1,3} = 2 Using Progression 1; // Vertical Lines
Transfinite Line {2,4} = num_cells Using Progression 1; // Horizontal lines
Transfinite Surface {1};
Recombine Surface {1};

// Define physical groups for boundary conditions
Physical Line("inlet", 1) = {1};
Physical Line("walls", 2) = {2, 4};
Physical Line("outlet", 3) = {3};

Physical Surface("fluid", 100) = {1};

Mesh 2;
Save "case.msh";
