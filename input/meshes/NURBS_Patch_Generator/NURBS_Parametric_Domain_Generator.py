
"""
**************************************************************
			NURBS Parametric Domain Generator (2D)
**************************************************************

This python script will be used to create a single NURBS patch to be used as a 
parametric domain in the DPGSolver. The knot domain (xi, eta) must be in [-1,1] 
(this is because the parametric domain in the DPGSolver code is in this domain). 

The setup process in the code will involve reading the knot vectors (xi and 
eta) as well as the design control point locations. The DPG code will then 
read in a GMSH file consisting of a reference domain (a square in the domain [-1,1] 
in both dimensions). This reference square will be discretized into the elements (on
the parametric domain). Using the NURBS mapping, each element's vertices and face nodes
will be found on the physical domain using their corresponding values on the 
parametric domain.

NOTE: 
	The output file will contain a section with all the control points listed as well
	as a connectivity section. This is done because it is possible for a patch to use
	the same control points (if the spline is closed for instance)

Optimization:
	- The NURBS patch will be able to be used for optimization as well. To ease the process
		and avoid the use of RBF, the patch will be such that there will be no RBF 
		volume points (no interior NURBS control points in the domain). The elements will
		not cross since the NURBS control points will be spaced apart and the parametric 
		domain will ensure the elements do not attain negative volumes.
"""

import Basis

import User_Defined_Patch
import Airfoil_Patch
import Internal_Channel_Patch

import numpy
import matplotlib.pyplot as plt
import sys

# The type of patch to use. Perhaps take this as a command line
# argument
#CONST_Patch_Type = "Internal_Channel_Patch"
#CONST_Patch_Type = "User_Defined_Patch"
CONST_Patch_Type = "Airfoil_Patch"
CONST_EPS = 1E-9

CONST_Output_file_name = "geometry_parameters_airfoil_P3_NumPtsXi12_Q1_NumPtsEta2.geo"

def NURBS_patch(xi,eta,BasisFunctionsList, ControlPoints_and_Weights, grad_index=None):

	"""
	Knowing the basis functions and locations of the points,
	compute the parametric NURBS patch expression

	:param xi: The xi value to evaluate the patch at
	:param eta: The eta value to evaluate the patch at
	:param BasisFunctionsList: The list (list of lists in 2D) of basis functions
	:param ControlPoints_and_Weights: The list (list of lists in 2D) of Points. 
		The i,j index of the Point structure points to a list of the form [x,y,weight], 
		which correspond to the physical location and weight of a given control point
	:param grad_index: An index value of 0 or 1 if the basis functions list provided 
		contains a matrix of lists holding the gradients of the basis functions

	:return: Value of the patch on the physical domain at the given xi,eta point on the 
		parametric domain. Using this function in a lambda expression will 
		ease the plotting process. The return list has two dimensions [x, y]. If the 
		gradients are provided, then the return list is [del_x, del_y], where del is the partial
		with respect to whichever parameter (xi or eta) is specified by grad_index
	"""

	physical_val = [0,0]

	numI = len(BasisFunctionsList)
	numJ = len(BasisFunctionsList[0])

	for i in range(numI):
		for j in range(numJ):
			for k in range(2):

				if grad_index is not None:
					physical_val[k] = physical_val[k] + BasisFunctionsList[i][j][grad_index](xi,eta) * ControlPoints_and_Weights[i][j][k]

				else:
					physical_val[k] = physical_val[k] + BasisFunctionsList[i][j](xi,eta) * ControlPoints_and_Weights[i][j][k]

	return physical_val


def plot_patch(patch_parameters):

	"""
	Plot the given patch using matplotlib. Display the control points as well
	as the constant knot lines. For the xi and eta
	values on the edge (boundary), subtract or add a small epsilon value to ensure
	they are in the patch domain.

	:param patch_parameters: The dictionary with the patch_parameters
	"""

	xiVector = patch_parameters["xiVector"]
	etaVector = patch_parameters["etaVector"]
	ControlPoints_and_Weights = patch_parameters["ControlPoints_and_Weights"]
	P = patch_parameters["P"]
	Q = patch_parameters["Q"]

	# Get the NURBS basis functions associated with each control point
	NURBS_basis_functions = Basis.get_NURBS_basis_functions(ControlPoints_and_Weights, P, Q, xiVector, etaVector)

	# Get the patch parametric function
	NURBS_patch_function = lambda xi,eta,BasisFunctionsList=NURBS_basis_functions, \
			ControlPoints_and_Weights=ControlPoints_and_Weights: \
			NURBS_patch(xi,eta, BasisFunctionsList, ControlPoints_and_Weights)

	# The distinct knot values along both parametric coordinates (xi and eta)
	xi_distinct_values = []
	eta_distinct_values = []

	for xi_val in xiVector:
		if xi_val not in xi_distinct_values:
			xi_distinct_values.append(xi_val)

	for eta_val in etaVector:
		if eta_val not in eta_distinct_values:
			eta_distinct_values.append(eta_val)

	# Discretize the xi and eta domains so we can plot the lines
	num_linspace_pts = 40
	eta_linspace = numpy.linspace(etaVector[0]+CONST_EPS, etaVector[-1]-CONST_EPS, num_linspace_pts)
	xi_linspace = numpy.linspace(xiVector[0]+CONST_EPS, xiVector[-1]-CONST_EPS, num_linspace_pts)

	# Plot the constant xi values:
	for xi_val in xi_distinct_values:
		if xi_distinct_values.index(xi_val) == 0:
			xi_val += CONST_EPS
		elif xi_distinct_values.index(xi_val) == (len(xi_distinct_values)-1):
			xi_val -= CONST_EPS

		xPlotPts = []
		yPlotPts = []

		for eta_val in eta_linspace:
			xy_pt = NURBS_patch_function(xi_val, eta_val)
			xPlotPts.append(xy_pt[0])
			yPlotPts.append(xy_pt[1])

		plt.plot(xPlotPts, yPlotPts, c='b')

	# Plot the constant eta values
	for eta_val in eta_distinct_values:
		if eta_distinct_values.index(eta_val) == 0:
			eta_val += CONST_EPS
		elif eta_distinct_values.index(eta_val) == (len(eta_distinct_values)-1):
			eta_val -= CONST_EPS

		xPlotPts = []
		yPlotPts = []

		for xi_val in xi_linspace:
			xy_pt = NURBS_patch_function(xi_val, eta_val)
			xPlotPts.append(xy_pt[0])
			yPlotPts.append(xy_pt[1])

		plt.plot(xPlotPts, yPlotPts, c='b')

	# Scatter the control points
	xPlotPts = []
	yPlotPts = []
	for i in range(len(ControlPoints_and_Weights)):
		for j in range(len(ControlPoints_and_Weights[0])):
			xPlotPts.append(ControlPoints_and_Weights[i][j][0])
			yPlotPts.append(ControlPoints_and_Weights[i][j][1])

	plt.scatter(xPlotPts, yPlotPts, c='r')

	plt.grid()
	plt.show(block=True)


def output_file(patch_parameters, output_file_name):

	"""
	Output the file with the patch parameters. The DPG code will
	read this file and use the values to perform the parametric
	mapping

	NOTE: Optimization_ControlPoints_and_Weights is a key that, when present, 
		corresponds to the control points used for the optimization. If this 
		data is present, then print it as well in the .geo file

	:param patch_parameters: The parameters used to define the patch
	"""

	xiVector = patch_parameters["xiVector"]
	etaVector = patch_parameters["etaVector"]
	ControlPoints_and_Weights = patch_parameters["ControlPoints_and_Weights"]
	P = patch_parameters["P"]
	Q = patch_parameters["Q"]

	num_xi_pts = len(ControlPoints_and_Weights)
	num_eta_pts = len(ControlPoints_and_Weights[0])

	# Load all the control points into a list (a connecitivity format will
	# be used to load the points)
	ControlPoints_and_Weights_list = []
	for i in range(num_xi_pts):
		for j in range(num_eta_pts):
			if ControlPoints_and_Weights[i][j] not in ControlPoints_and_Weights_list:
				ControlPoints_and_Weights_list.append(ControlPoints_and_Weights[i][j])

	with open(output_file_name, "w") as fp:

		fp.write("/** Geometry parameters for test case: euler/steady/NURBS\n")
		fp.write("*/\n\n")

		# The order of the patch in each parameter direction
		fp.write("P(xi_order) %d \n" % P)
		fp.write("Q(eta_order) %d \n" % Q)
		fp.write("\n")

		# The knot vectors
		fp.write("knots_xi %d \n" % len(xiVector))
		for val in xiVector:
			fp.write("%.14e \n" % val)
		fp.write("\n")
		fp.write("knots_eta %d \n" % len(etaVector))
		for val in etaVector:
			fp.write("%.14e \n" % val)
		fp.write("\n")
		fp.write("\n")

		# Control Point Data information
		fp.write("Control_Point_Data %d \n" % len(ControlPoints_and_Weights_list))
		for pt in ControlPoints_and_Weights_list:
			fp.write("%.14e %.14e %.14e \n" % (pt[0], pt[1], pt[2]))
		fp.write("\n")

		# Connectivity Information	
		fp.write("Control_Point_Connectivity %d %d\n" % (num_xi_pts, num_eta_pts))
		for i in range(num_xi_pts):
			for j in range(num_eta_pts):
				fp.write("%d " % ControlPoints_and_Weights_list.index(ControlPoints_and_Weights[i][j]))
			fp.write("\n")
		fp.write("\n")

		if "Optimization_ControlPoints_and_Weights" in patch_parameters:
			
			#Optimization Information
			optimization_data_tuples = patch_parameters['Optimization_ControlPoints_and_Weights']

			fp.write("Optimization_Point_Connectivity %d\n" % (len(optimization_data_tuples)))
			for data_tuple in optimization_data_tuples:
				fp.write("%d %d %d \n" % (ControlPoints_and_Weights_list.index(data_tuple[0]), data_tuple[1], data_tuple[2]))
			fp.write("\n")	

			optimization_control_pt_limit_tuples = patch_parameters['Optimization_ControlPoints_Limits']
			fp.write("Optimization_Point_Limit %d\n" % (len(optimization_control_pt_limit_tuples)))
			for data_tuple in optimization_control_pt_limit_tuples:
				fp.write("%d %e %e \n" % (ControlPoints_and_Weights_list.index(data_tuple[0]), data_tuple[1], data_tuple[2]))
			fp.write("\n")	


		if "area_ref" in patch_parameters:
			fp.write("area_ref = %e;" % patch_parameters['area_ref'])
			fp.write("\n")

		if "cm_le_x" in patch_parameters:
			fp.write("cm_le_x = %e;" % patch_parameters['cm_le_x'])
			fp.write("\n")

		if "cm_le_y" in patch_parameters:
			fp.write("cm_le_y = %e;" % patch_parameters['cm_le_y'])
			fp.write("\n")


def test():

	"""
	Test the patch and some parameters in it
	"""

	# Get the patch parameters
	patch_parameters = user_defined_patch.get_patch_information()

	xiVector = patch_parameters["xiVector"]
	etaVector = patch_parameters["etaVector"]
	ControlPoints_and_Weights = patch_parameters["ControlPoints_and_Weights"]
	P = patch_parameters["P"]
	Q = patch_parameters["Q"]

	# Get the NURBS basis functions associated with each control point
	NURBS_basis_functions = Basis.get_NURBS_basis_functions(ControlPoints_and_Weights, P, Q, xiVector, etaVector)

	# Get the patch parametric function
	NURBS_patch_function = lambda xi,eta,BasisFunctionsList=NURBS_basis_functions, \
			ControlPoints_and_Weights=ControlPoints_and_Weights: \
			NURBS_patch(xi,eta, BasisFunctionsList, ControlPoints_and_Weights)

	
	# Get the patch gradients
	NURBS_basis_functions_grad = Basis.get_grad_NURBS_basis_functions(ControlPoints_and_Weights, P, Q, xiVector, etaVector)

	NURBS_patch_function_del_xi = lambda xi,eta,BasisFunctionsList=NURBS_basis_functions_grad, \
			ControlPoints_and_Weights=ControlPoints_and_Weights: \
			NURBS_patch(xi,eta, BasisFunctionsList, ControlPoints_and_Weights, 0)

	NURBS_patch_function_del_eta = lambda xi,eta,BasisFunctionsList=NURBS_basis_functions_grad, \
			ControlPoints_and_Weights=ControlPoints_and_Weights: \
			NURBS_patch(xi,eta, BasisFunctionsList, ControlPoints_and_Weights, 1)
	
	xi_test = -0.25
	eta_test = 0.5

	print("C(xi, eta)    : " + str(NURBS_patch_function(xi_test, eta_test)))
	print("C_xi(xi,eta)  : " + str(NURBS_patch_function_del_xi(xi_test, eta_test)))
	print("C_eta(xi,eta) : " + str(NURBS_patch_function_del_eta(xi_test, eta_test)))
	
	return

	# Finite difference tests
	h = 1E-8
	
	del_x_del_xi_fd = (1./h) * (NURBS_patch_function(xi_test + h, eta_test)[0] - NURBS_patch_function(xi_test, eta_test)[0])
	del_y_del_xi_fd = (1./h) * (NURBS_patch_function(xi_test + h, eta_test)[1] - NURBS_patch_function(xi_test, eta_test)[1])

	del_x_del_eta_fd = (1./h) * (NURBS_patch_function(xi_test, eta_test + h)[0] - NURBS_patch_function(xi_test, eta_test)[0])
	del_y_del_eta_fd = (1./h) * (NURBS_patch_function(xi_test, eta_test + h)[1] - NURBS_patch_function(xi_test, eta_test)[1])
	
	print("C_xi_fd(xi,eta)  : " + str([del_x_del_xi_fd, del_y_del_xi_fd]))
	print("C_eta_fd(xi,eta) : " + str([del_x_del_eta_fd, del_y_del_eta_fd]))
	

def main():
	
	"""
	Main method of the patch generation module
		
	Command line argument options:
		plot = Plot the given patch to see what it looks like. This will
			plot the outer lines of the patch (-1, 1 constant xi and eta lines 
			for the 2D case)
		output_file = Output the parametric NURBS patch file which will be able
			to be processed by the DPG code.
		[No Arguments] = Run the plot and output_file case (together)

	:return : -

	"""

	# Get the patch parameters
	if CONST_Patch_Type == "User_Defined_Patch":
		patch_parameters = User_Defined_Patch.get_patch_information()
	elif CONST_Patch_Type == "Airfoil_Patch":
		patch_parameters = Airfoil_Patch.get_patch_information(P=3, num_ctrl_pts_xi=12)
	elif CONST_Patch_Type == "Internal_Channel_Patch":
		patch_parameters = Internal_Channel_Patch.get_patch_information()
	else:
		raise ValueError("Unknown Patch Type")
	output_file_name = "geometry_parameters_airfoil_P3_NumPtsXi12_Q1_NumPtsEta2.geo"

	# Parse command line arguments
	if len(sys.argv) == 1:
		# No command line arguments so plot and output the file
		plot_patch(patch_parameters)
		output_file(patch_parameters, output_file_name)

	elif(sys.argv[1] == "plot"):
		# Only plot the patch
		plot_patch(patch_parameters)
	
	elif(sys.argv[1] == "output_file"):
		# Only output the patch file
		output_file(patch_parameters, output_file_name)


if __name__ == "__main__":
	main()
	#test()




