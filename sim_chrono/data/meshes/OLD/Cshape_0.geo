SetFactory("OpenCASCADE");

w = DefineNumber[ 0.02, Name "Parameters/w" ];
//+
h = DefineNumber[ 0.02, Name "Parameters/h" ];
//+
fillet_r = DefineNumber[ 0.005, Name "Parameters/fillet_r" ];
//+
cs_rad = DefineNumber[ 0.001, Name "Parameters/cs_rad" ];
//+
gap_h = DefineNumber[ 0.003, Name "Parameters/gap_h" ];
//+
gap_angle = DefineNumber[ 0.0, Name "Parameters/gap_angle" ];
//+
Point(1) = {w, gap_h, 0, 1.0};
//+
Point(2) = {w-fillet_r, h-fillet_r, 0, 1.0}; // center
//+
Point(3) = {w, h-fillet_r, 0, 1.0}; // first point
//+
Point(4) = {w-fillet_r, h, 0, 1.0}; // second point
//+
Point(5) = {fillet_r, h-fillet_r, 0, 1.0}; // center
//+
Point(6) = {fillet_r, h, 0, 1.0}; // first point
//+
Point(7) = {0.0, h-fillet_r, 0, 1.0}; // second point
//+
Point(8) = {0.0, 0.0, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {4, 6};
//+
Line(3) = {7, 8};
//+
Circle(6) = {0.0, 0.0, 0, cs_rad, 0, 2*Pi};
//+
Rotate {{1, 0, 0}, {0.0, 0.0, 0}, Pi/2} {
  Curve{6}; 
}
//+
Circle(7) = {3, 2, 4};
//+
Circle(8) = {6, 5, 7};
//+
Extrude {{0, 1, 0}, {0, 0, 0}, Pi/4} {
  Curve{6}; Layers{5}; Recombine;
}
//+
Extrude {0, 0.01, 0} {
  Curve{6}; Layers{5}; Recombine;
}
