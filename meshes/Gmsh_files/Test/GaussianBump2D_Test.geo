// Gmsh project created on Thu Aug 31 16:19:00 2017
//+
Point(1) = {-4, 0, 0, 1.0};
//+
Point(2) = {4, 0, 0, 1.0};
//+
Point(3) = {4, 1.5, 0, 1.0};
//+
Point(4) = {-4, 1.5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
