Include "../parameters.geo";
//mesh_level = 2; mesh_type = MIXED; mesh_domain = STRAIGHT; pde_name = ADVECTION; geom_adv = GEOM_ADV_XL; geom_unaligned = 1;
//mesh_level = 2; mesh_type = QUAD; mesh_domain = PARAMETRIC; pde_name = EULER; geom_adv = GEOM_ADV_INTERNAL; geom_unaligned = 0; geom_ar = 0.5;

// gmsh file for the NURBS patch airfoil case (creates a stretched grid)

// Geometry Specification
l = 1;
h = 1;

Point(1) = {-l,-h,-0,lc};
Point(2) = {+l,-h,-0,lc};
Point(3) = {-l,+h,-0,lc};
Point(4) = {+l,+h,-0,lc};
Point(5) = {-0,-h,-0,lc};
Point(6) = {+0,+h,-0,lc};

Line(1001) = {1,5};
Line(1002) = {5,2};
Line(1003) = {3,6};
Line(1004) = {6,4};
Line(2001) = {1,3};
Line(2002) = {2,4};
Line(2003) = {5,6};


If (mesh_level == 0)
	
	Transfinite Line{1001:1002} = 4+1 Using Progression 1;
	Transfinite Line{2001:2002} = 5+1 Using Progression 1.1;
	Transfinite Line{2003}      = 5+1 Using Progression 1.1;
	Transfinite Line{1003:1004} = 4+1 Using Progression 1;

ElseIf (mesh_level == 1)

	Transfinite Line{1001:1002} = PYTH_TEMP_NUM_Y_ELEMENTS Using Progression 1;
	Transfinite Line{2001:2002} = PYTH_TEMP_NUM_X_ELEMENTS Using Progression 1.15;
	Transfinite Line{2003}      = PYTH_TEMP_NUM_X_ELEMENTS Using Progression 1.15;
	Transfinite Line{1003:1004} = PYTH_TEMP_NUM_Y_ELEMENTS Using Progression 1;

ElseIf (mesh_level == 2)

	Transfinite Line{1001:1002} = 16+1 Using Progression 1;
	Transfinite Line{2001:2002} = 20+1 Using Progression 1.2;
	Transfinite Line{2003}      = 20+1 Using Progression 1.2;
	Transfinite Line{1003:1004} = 16+1 Using Progression 1;

ElseIf (mesh_level == 3)

	Transfinite Line{1001:1002} = 44+1 Using Progression 1;
	Transfinite Line{2001:2002} = 55+1 Using Progression 1.15;
	Transfinite Line{2003}      = 55+1 Using Progression 1.15;
	Transfinite Line{1003:1004} = 44+1 Using Progression 1;

ElseIf (mesh_level == 4)

	Transfinite Line{1001:1002} = 64+1 Using Progression 1;
	Transfinite Line{2001:2002} = 80+1 Using Progression 1.15;
	Transfinite Line{2003}      = 80+1 Using Progression 1.15;
	Transfinite Line{1003:1004} = 64+1 Using Progression 1;

Else
	Transfinite Line{1001:1002} = 2^(mesh_level)*2+1 Using Progression 1;
	Transfinite Line{2001:2002} = 2^(mesh_level+1)+1 Using Progression 1.3;
	Transfinite Line{2003}      = 2^(mesh_level+1)+1 Using Progression 1.3;
	Transfinite Line{1003:1004} = 2^(mesh_level)*2+1 Using Progression 1;
EndIf

Line Loop (4001) = {1001,2003,-1003,-2001};
Line Loop (4002) = {1002,2002,-1004,-2003};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface {4001} AlternateLeft;
Transfinite Surface {4002} AlternateLeft;

If (mesh_type == TRI)
	// Do nothing.
ElseIf (mesh_type == QUAD)
	Recombine Surface{4001,4002};
ElseIf (mesh_type == MIXED)
	Recombine Surface{4002};
Else
	Error("Unsupported mesh_type: %d",mesh_type); Exit;
EndIf

// Physical parameters for '.msh' file
bc_straight =   BC_STEP_SC;
bc_curved   = 2*BC_STEP_SC;
If (mesh_domain == STRAIGHT)
	bc_base = bc_straight;
Else
	bc_base = bc_curved;
EndIf

If (pde_name == EULER)
	If (geom_adv == GEOM_NURBS_AIRFOIL_O_GRID)
		// Periodic left and right x face
		//Physical Line(bc_base+BC_BACKPRESSURE) = {2001};
		//Physical Line(bc_base+BC_TOTAL_TP) = {2002};
		
		Physical Line(bc_base+PERIODIC_XL) = {2001};
	 	Physical Line(bc_base+PERIODIC_XR) = {2002};
		Periodic Line {2002} = {2001}; // Periodic (x)

		// Slipwall for bottom face and reimann for top
		Physical Line (bc_base+BC_SLIPWALL)= {1001:1002};
		Physical Line (bc_base+BC_RIEMANN) = {1003:1004};
	Else
		Error("Unsupported geom_adv: %d",geom_adv); Exit;
	EndIf
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

Physical Surface(9401) = 4001;
Physical Surface(9402) = 4002;


// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
