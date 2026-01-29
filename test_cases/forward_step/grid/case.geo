// This GMSH file is drawn in units of meters

boundary_layer_height = 0.025; // Height of the boundary layer region
step_height = 0.2; // Height of the forward step

// Define points starting in lower left corner and going
// clockwise around the forward-step geometry
Point(1) = {0.0, 0.0, 0.0}; // Bottom left corner
Point(2) = {0.0, 1.0, 0.0}; // Top left corner
Point(3) = {3.0, 1.0, 0.0}; // Top right corner
Point(4) = {3.0, step_height+boundary_layer_height, 0.0}; // Boundary layer height
Point(5) = {3.0, step_height, 0.0};
Point(6) = {0.6, step_height, 0.0};
Point(7) = {0.6, 0.0, 0.0};

Point(9) = {0.6, step_height+boundary_layer_height, 0.0}; // Point above the step corner

// Define lines that make up the domain
Line(1) = {1, 2}; // Inlet
Line(2) = {2, 3}; // Top wall
Line(3) = {3, 4}; // Outlet above boundary layer
Line(4) = {4, 5}; // Outlet inside the boundary layer
Line(5) = {5, 6}; // Top of step
Line(6) = {6, 7}; // Front of step
Line(7) = {7, 1}; // Bottom wall upstream of the step
Line(9) = {6, 9}; // Interface: vertical edge of boundary-layer block
Line(10) = {9, 4}; // Interface: horizontal edge of boundary-layer block


// Only the boundary layer region is structured/transfinite
Curve Loop(1) = {9, 10, 4, 5}; // Upstream, upper section
Surface(1) = {-1};
Transfinite Line {9,4} = 80 Using Progression 1; // Vertical Lines
Transfinite Line {10,5} = 200 Using Progression 1; // Horizontal lines
Transfinite Surface {1};
Recombine Surface {1};

// Rest of the domain is unstructured (shares interface with Surface 1)
Curve Loop(2) = {1, 2, 3, -10, -9, 6, 7}; // Outer domain
Plane Surface(2) = {2};

Transfinite Line {1} = 150 Using Progression 1; // Inlet
Transfinite Line {2} = 200 Using Progression 1; // Top wall
Transfinite Line {3} = 50 Using Progression 1; // Outlet above boundary layer
Transfinite Line {6} = 60 Using Progression 1; // Front of step
Transfinite Line {7} = 60 Using Progression 1; // Bottom wall upstream of the step


Physical Line("inlet", 1) = {1};
Physical Line("walls", 2) = {2, 5, 6, 7};
Physical Line("outlet", 3) = {3, 4};

Physical Surface("fluid", 100) = {1, 2};

Mesh 2;
Save "case.msh";
