// This GMSH file is drawn in millimeters.
// The bump has a chord length of 1000mm and the height of the bump is 4.2% of the bump chord length.

height_factor = 0.042; // Height of bump as a fraction of chord length
chord = 1000; // Chord length of bump (mm)
bump_height = height_factor*chord; // Height of bump


// Center of the bump
CX = 0.0; // X coordinate of leading edge of bump
CY = 0.0; // Y coordinate of leading edge of bump

// Outer domain parameters
upstream_spacing = 1000; // Distance from the bump to the inlet
downstream_spacing = 2000; // Distance from the bump to the outlet
upper_spacing = 1000; // Distance from the bump to the top wall


// Points for the bump
Point(1) = {CX, CY, 0}; // Leading edge of bump
Point(2) = {CX + 0.5*chord, CY+bump_height, 0}; // Upper point of bump
Point(3) = {CX + chord, CY, 0}; // Trailing edge of bump

// This is the center of a circle that creates a 30 degree angle with respect to the vertical
Radius = (chord*chord + 4*bump_height*bump_height) / (8*bump_height); // Radius of the circle
Printf("Radius = %f", Radius);
Point(4) = {CX + 0.5*chord, CY-(Radius-bump_height), 0}; // Center point used to create the bump arc

// Points for rectangular outer domain
Point(5) = {CX - upstream_spacing, CY, 0}; // Upstream lower point
Point(6) = {CX - upstream_spacing, CY + upper_spacing, 0};   // Upstream upper point
Point(7) = {CX + downstream_spacing, CY + upper_spacing, 0}; // Downstream upper point
Point(8) = {CX + downstream_spacing, CY, 0}; // Downstream lower point

Circle(1) = {1, 4, 2}; // Upstream arc of bump
Circle(2) = {3, 4, 2}; // Downstream arc of bump

// Lines connecting outer domain
Line(3) = {1, 5}; // Bottom upstream
Line(4) = {5, 6}; // Upstream vertical line (inlet)
Line(5) = {6, 7}; // Top wall
Line(6) = {7, 8}; // Downstream vertical line (outlet)
Line(7) = {8, 3}; // Bottom downstream


// Surface
Curve Loop(1) = {3, 4, 5, 6, 7, 2, -1};
Plane Surface(1) = {1};


//---- Sizing field for uniform unstructured mesh near bump----
// Create a Distance field
DistanceField = 1;
Field[DistanceField] = Distance;
Field[DistanceField].CurvesList = {1, 2};
Field[DistanceField].Sampling = 100;

// Apply a Threshold field to control mesh size growth
ThresholdField = 2;
Field[ThresholdField] = Threshold;
Field[ThresholdField].IField = DistanceField; // Use the Distance field
Field[ThresholdField].LcMin = 5; // Minimum cell size
Field[ThresholdField].LcMax = 35; // Maximum cell size
Field[ThresholdField].DistMin = 0; // Distance for min cell size
Field[ThresholdField].DistMax = 1000; // Distance for max cell size

// Set the Background Field
Background Field = ThresholdField;
//-------------------------------------------------------------------


// These numbers are the boundary condition tags
Physical Line("inlet", 1) = {4};
Physical Line("bump", 2) = {1, 2};
Physical Line("top_wall", 3) = {5};
Physical Line("outlet", 4) = {6};
Physical Line("symmetry", 5) = {3,7};

// Physical surface tag must be unique vs line tags
Physical Surface("fluid", 100) = {1};

Mesh 2;
Save "case.msh";
