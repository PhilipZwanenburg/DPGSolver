"""
Module: Airfoil_Patch.py
------------------------------------------

Creates a NURBS patch around an airfoil. Use a least squares approach to determine the 
locations of the control points. Take as input the number of control points desired and 
spline order. 

Use a standard approach to place the slipwall on one xi/eta location, 
farfield on another, etc...

"""

import math
import numpy
import matplotlib.pyplot as plt
import Basis
import sys
import scipy.integrate


# ==================================================
# 				Patch Parameters

# The properties of the patch in the xi direction (which
# traverses around the airfoil from the trailing edge, bottom
# surface to the leading edge and back)
CONST_P = 2
CONST_NUM_CONTROL_PTS_XI = 25

# Properties of the patch in the eta direction (eta increases in the 
# normal direction from the airfoil surface to the farfield)
CONST_Q = 1
CONST_ETA_KNOTS = [-1, -1, 1, 1]


CONST_R_FARFIELD = 20.

# ==================================================


CONST_Airfoil_Type = "NACA0012"
#CONST_Airfoil_Type = "NACA4412"

CONST_QUADRATURE_P = 35

# Plot Parameters (for testing)
CONST_PlotXRange = [-0.7, 0.7]
CONST_PlotYRange = [-0.3, 0.3]


def airfoil_parametric_equation(s):

	"""
	Get the coordinates of the airfoil (x,y) based on the 
	parametric coordinate s. s is an element of [-1, 1], and will
	go from the trailing edge, through the lower surface of the airfoil
	to the leading edge (s = 0) and back over the upper surface.
	The x limit of the airfoil is from -0.5 to 0.5. 

	:param s: The parametric coordinate

	:return: A tuple for the x and y coordinate for the point on the
		airfoil
	"""

	if s >= -1 and s <= 0:
		x_by_c = -1*s

		x = x_by_c*1.0 + -0.5  # Chord length  = 1.0
		y = NACA_equation(x_by_c)[1]

	else:
		x_by_c = s

		x = x_by_c*1.0 + -0.5  # Chord length  = 1.0
		y = NACA_equation(x_by_c)[0]

	return (x,y)


def NACA_equation(x_by_c):

	"""
	Get the top and bottom y value for the given airfoil by
	using the NACA Airfoil equations. 

		m = The max value of the mean line in hundreths of chord
		p = The chordwise position of the maximum chamber in tenths
			of chord
		t = Max thickness normalized by the chord length

	:param x_by_c: The normalized location (by chord length) to
		evaluate properties at

	:return : Tuple containing the top and bottom location of the 
		airfoil surface at the given chord length point. 
		(top_y_value, bottom_y_value)
	"""

	if CONST_Airfoil_Type == "NACA0012":
		m = 0.0
		p = 0.0
		t = 0.12
	elif CONST_Airfoil_Type == "NACA4412":
		m = 0.04
		p = 0.4
		t = 0.12
	else:
		raise ValueError("Not Supported Airfoil")

	x = x_by_c

	# Calculate y_c
	y_c = 0
	if not (p == 0 and m == 0):
		if x < p:
			y_c = (m/(p**2.))*(2.*p*x - x**2.)
		else:
			y_c = (m/((1.-p)**2.))*(1. - 2.*p + 2.*p*x - x**2.)

	# Calculate y_t

	a0 = 1.4845
	a1 = -0.6300
	a2 = -1.7580
	a3 = 1.4215
	a4 = -0.5180

	y_t = (t/1.0)*(a0*math.sqrt(x) + a1*x + a2*x**2. + a3*x**3. + a4*x**4.)

	return y_c + y_t, y_c - y_t


def integrate_gaussian_quadrature(function, nodes, weights):

	"""
	Integrate the function provided on the domain [-1,1] using 
	gaussian quadrature

	:param function: The lambda expression for the function
	:param nodes: The nodes for the quadrature
	:param weights: The weights for the quadrature

	:return: The value of the integral of the given function from -1 to 1
	"""

	val = 0
	for i in range(len(nodes)):
		val += function(nodes[i]) * weights[i]

	return val


def get_BSpline_control_points_continuous_least_square(BSpline_Basis, knots, n, ref_function):

	"""
	Use a continuous least square approximation to get the location of the 
	control points. 

	NOTE: An open knot vector is always used for the approximation

	:param BSpline_Basis: The list of lambda expressions for the B Spline basis functions
	:param knots: The knot vector
	:param n: The number of control points to find
	:param ref_function: The parametric curve to approximate (Must have a domain of [-1,1])
	
	:return: The the list of control points for the B Spline
	"""
	
	# Create the list to hold the control points. Since an open knot vector
	# is used, the first and last control points are identical to the first and 
	# last airfoil points
	Control_Pts = []
	for i in range(n):
		if i == 0:
			Control_Pts.append(ref_function(-1.))
		elif i == n-1:
			Control_Pts.append(ref_function(1.))
		else:
			Control_Pts.append([None, None])


	# The nodes and weights for the quadrature
	nodes, weights = numpy.polynomial.legendre.leggauss(CONST_QUADRATURE_P)

	# Least square is of the form [M]*{X} = {b}. [M] matrix is the 
	# same for the x and y least square computation

	# Form the M matrix
	M = numpy.zeros((n-2, n-2))

	for i in range(n-2):
		
		# Basis functions range from i=2 to i=(n-1) (starting from i=1)
		N_ip = BSpline_Basis[i+1]  

		for j in range(n-2):
			
			# Basis functions range from j=2 to j=(n-1) (starting from j=1)
			N_jp = BSpline_Basis[j+1]
			
			func = lambda xi, N_ip=N_ip, N_jp=N_jp: N_ip(xi)*N_jp(xi)
			M[i][j] = integrate_gaussian_quadrature(func, nodes, weights)

	# Form the b vector (for x and y coordinate least square)
	b_x = numpy.zeros((n-2, 1))
	b_y = numpy.zeros((n-2, 1))

	for j in range(n-2):
		
		# Basis functions range from j=2 to j=(n-1) (starting from j=1)
		N_jp = BSpline_Basis[j+1]

		N_1p = BSpline_Basis[0]
		N_np = BSpline_Basis[n-1]

		P1x = Control_Pts[0][0]
		P1y = Control_Pts[0][1]

		Pnx = Control_Pts[n-1][0]
		Pny = Control_Pts[n-1][1]

		term1_x = integrate_gaussian_quadrature(
			lambda xi, N_jp=N_jp: ref_function(xi)[0]*N_jp(xi), 
			nodes, 
			weights)
		
		term1_y = integrate_gaussian_quadrature(
			lambda xi, N_jp=N_jp: ref_function(xi)[1]*N_jp(xi), 
			nodes, 
			weights)

		term2 = integrate_gaussian_quadrature(
			lambda xi, N_jp=N_jp, N_1p=N_1p: N_1p(xi)*N_jp(xi), 
			nodes, 
			weights)

		term3 = integrate_gaussian_quadrature(
			lambda xi, N_jp=N_jp, N_np=N_np: N_np(xi)*N_jp(xi), 
			nodes, 
			weights)

		b_x[j][0] = term1_x - P1x * term2 - Pnx * term3
		b_y[j][0] = term1_y - P1y * term2 - Pny * term3

	control_pts_x = numpy.linalg.solve(M, b_x)
	control_pts_y = numpy.linalg.solve(M, b_y)

	# Load the computed control pts into the control points array.
	for i in range(n-2):
		Control_Pts[i+1] = (control_pts_x[i][0], control_pts_y[i][0])

	return Control_Pts


def get_BSpline_control_points_discrete_least_square(BSpline_Basis, knots, n):

	"""
	Use a least square approximation to get the location of the control points.
	In this method, a discrete least square approximation will be done. So, a set
	of points will be found on the airfoil and then these will be used to get
	the B Spline control points.

	NOTE: An open knot vector is always used for the approximation

	:param BSpline_Basis: The list of lambda expressions for the B Spline basis functions
	:param knots: The knot vector
	:param n: The number of control points to find
	
	:return: The the list of control points for the B Spline
	"""

	m = n+10  # The number of airfoil points to use for the least square
	t_min = knots[0]
	t_max = knots[-1]
	t_k = numpy.linspace(t_min, t_max, m)

	# Get the list of airfoil points (used as the points to approximate in a 
	# least square sense)
	Q = []
	for t in t_k:
		Q.append(airfoil_parametric_equation(t))

	# Create the list to hold the control points. Since an open knot vector
	# is used, the first and last control points are identical to the first and 
	# last airfoil points
	Control_Pts = []
	for i in range(n):
		if i == 0:
			Control_Pts.append(Q[0])
		elif i == n-1:
			Control_Pts.append(Q[-1])
		else:
			Control_Pts.append([None, None])


	# Least square is of the form [M]*{X} = {b}. [M] matrix is the 
	# same for the x and y least square computation

	# Form the M matrix
	M = numpy.zeros((n-2, n-2))

	for i in range(n-2):
		
		# Basis functions range from i=2 to i=(n-1) (starting from i=1)
		N_ip = BSpline_Basis[i+1]  

		for j in range(n-2):
			
			# Basis functions range from j=2 to j=(n-1) (starting from j=1)
			N_jp = BSpline_Basis[j+1]
			term = 0

			for k in range(1, len(t_k)):
				
				# loop over the t_k values from t_2 to t_(m-1) (if t_k starts at t_1)
				term += N_ip(t_k[k])*N_jp(t_k[k])
			
			M[i][j] = term
	
	# Form the b vector (for x and y coordinate least square)
	b_x = numpy.zeros((n-2, 1))
	b_y = numpy.zeros((n-2, 1))

	for j in range(n-2):
		
		# Basis functions range from j=2 to j=(n-1) (starting from j=1)
		N_jp = BSpline_Basis[j+1]

		N_1p = BSpline_Basis[0]
		N_np = BSpline_Basis[n-1]

		P1x = Control_Pts[0][0]
		P1y = Control_Pts[0][1]

		Pnx = Control_Pts[n-1][0]
		Pny = Control_Pts[n-1][1]

		term1_x = 0
		term1_y = 0
		term2 = 0
		term3 = 0

		for k in range(1, len(t_k)):

			term1_x += Q[k][0] * N_jp(t_k[k])
			term1_y += Q[k][1] * N_jp(t_k[k])

			term2 += N_1p(t_k[k]) * N_jp(t_k[k])
			term3 += N_np(t_k[k]) * N_jp(t_k[k])

		b_x[j][0] = term1_x - P1x * term2 - Pnx * term3
		b_y[j][0] = term1_y - P1y * term2 - Pny * term3

	control_pts_x = numpy.linalg.solve(M, b_x)
	control_pts_y = numpy.linalg.solve(M, b_y)

	# Load the computed control pts into the control points array.
	for i in range(n-2):
		Control_Pts[i+1] = (control_pts_x[i][0], control_pts_y[i][0])

	sys.exit(0)

	return Control_Pts


def get_airfoil_BSpline_parameters():

	"""
	Get the parameters for the B Spline that will be used to 
	parametrize the airfoil surface. This method will generate the 
	knot vector and the control point locations using a least square
	approach.

	:return: Dictionary with the parameters for the B Spline
		return_dict = {
			'knots' : [list with the knots],
			'control_points': [[pt1_x, pt1_y], [pt2_x, pt2_y], ...],
			'P' : Order of the spline,
			'Spline_Function' : Lambda expression for the spline
		}
	"""

	# Spline parameters
	n = CONST_NUM_CONTROL_PTS_XI
	p = CONST_P
	m = n + p + 1  # number of elements in knot vector

	# Create a uniform open knot vector (use an open knot vector so first and last
	# knot are repeated p+1 times). Also, domain of the knot vector is [-1,1]
	knots = []
 
	for i in range(p+1):
		knots.append(-1.)

	num_non_end_knots = m - 2*(p+1)
	delta_knot = 2./(num_non_end_knots+1)

	for i in range(1, num_non_end_knots+1):
		knots.append(-1. + i*delta_knot)

	for i in range(p+1):
		knots.append(1.)


	# Get the B Spline basis functions
	BSpline_Basis = Basis.get_BSpline_basis_functions_1D(p, knots)

	# Compute the control points using a least square approximation
	#control_points = get_BSpline_control_points_discrete_least_square(BSpline_Basis, knots, n)
	control_points = get_BSpline_control_points_continuous_least_square(BSpline_Basis, knots, n, airfoil_parametric_equation)

	# Get the lambda expression for the spline
	spline_function = lambda xi, BSpline_Basis=BSpline_Basis, control_points=control_points: \
		Basis.get_Spline_Function(xi, control_points, BSpline_Basis)

	return {
		"knots" : knots,
		"control_points" : control_points,
		"P" : p,
		"Spline_Function" : spline_function
	}


def get_farfield_BSpline_parameters():

	"""
	Get the parameters for the B Spline that will be used to 
	parametrize the farfield surface. This method will generate the 
	knot vector and the control point locations using a least square
	approach. Since an o-grid is intended, this method will
	create a spline that approximates the farfield circle

	:return: Dictionary with the parameters for the B Spline
		return_dict = {
			'knots' : [list with the knots],
			'control_points': [[pt1_x, pt1_y], [pt2_x, pt2_y], ...],
			'P' : Order of the spline,
			'Spline_Function' : Lambda expression for the spline
		}
	"""

	# Spline parameters
	n = CONST_NUM_CONTROL_PTS_XI
	p = CONST_P
	m = n + p + 1  # number of elements in knot vector

	# Create a uniform open knot vector (use an open knot vector so first and last
	# knot are repeated q+1 times). Also, domain of the knot vector is [-1,1]
	knots = []
 
	for i in range(p+1):
		knots.append(-1.)

	num_non_end_knots = m - 2*(p+1)
	delta_knot = 2./(num_non_end_knots+1)

	for i in range(1, num_non_end_knots+1):
		knots.append(-1. + i*delta_knot)

	for i in range(p+1):
		knots.append(1.)

	# Create the lambda expression for the farfield (circle)
	farfield = lambda xi: (CONST_R_FARFIELD*math.cos(-1.*(xi+1)*math.pi), 
							CONST_R_FARFIELD*math.sin(-1.*(xi+1)*math.pi))

	# Get the B Spline basis functions
	BSpline_Basis = Basis.get_BSpline_basis_functions_1D(p, knots)

	# Compute the control points using a least square approximation
	control_points = get_BSpline_control_points_continuous_least_square(BSpline_Basis, knots, n, farfield)

	# Get the lambda expression for the spline
	spline_function = lambda xi, BSpline_Basis=BSpline_Basis, control_points=control_points: \
		Basis.get_Spline_Function(xi, control_points, BSpline_Basis)

	return {
		"knots" : knots,
		"control_points" : control_points,
		"P" : p,
		"Spline_Function" : spline_function
	}


def plot_airfoil_test():
	# Test method for plotting airfoil
	
	s = numpy.linspace(-1, 1, 201)

	x_vals = []
	y_vals = []

	for s_val in s:
		pt = airfoil_parametric_equation(s_val)
		x_vals.append(pt[0])
		y_vals.append(pt[1])

	plt.plot(x_vals, y_vals, c='b')


def plot_spline(spline_function, xi_range, control_points):

	xi_min = xi_range[0]
	xi_max = xi_range[-1]

	num_plot_pts = 101

	x_vals = []
	y_vals = []

	xi_vals = numpy.linspace(xi_min, xi_max, 201)

	for xi in xi_vals:
		pt = spline_function(xi)
		x_vals.append(pt[0])
		y_vals.append(pt[1])

	plt.plot(x_vals, y_vals, c='r')

	x_ctrl_pt = []
	y_ctrl_pt = []

	for pt in control_points:
		x_ctrl_pt.append(pt[0])
		y_ctrl_pt.append(pt[1])

	plt.scatter(x_ctrl_pt, y_ctrl_pt)


def get_optimization_pts(ControlPoints_and_Weights):

	"""
	Go through the Control Point net and set the points that will 
	act as the optimization control points. This method will combine
	the points into a list and also specify each point's degrees of 
	freedom. 

	:param ControlPoints_and_Weights: The net (matrix) of control points and
		weights.

	:return : A list with the control points for the optimization. The list
		will hold tuples of the form (cntrl_pt, dof_x_bool, dof_y_bool)
		where dof_x_bool = 0 if the x direction cannot be used as a degree of 
		freedom and is 1 if it can be (same for dof_y_bool)
	"""

	# For the Airfoil case, the trailing edge will be fixed and all other points
	# will be able to freely move in the y direction

	# The leading edge will also be fixed (although the spline does not interpolate
	# this location)

	optimization_control_pt_list = []

	num_airfoil_pts = len(ControlPoints_and_Weights) - 1  # Subtract one because first and last pt identical

	leading_edge_pt_index = int(len(ControlPoints_and_Weights)/2)

	for i in range(1,num_airfoil_pts):
		if i != leading_edge_pt_index:
			optimization_control_pt_list.append((ControlPoints_and_Weights[i][0], 0, 1))

	return optimization_control_pt_list


def get_patch_information():

	"""
	Return the Patch information

	:return : A dictionary holding key, value pairs
		for the patch information
		{
			xiVector: [],
			etaVector: [],
			ControlPoints_and_Weights: [],
			P : val,
			Q : val
		}
	"""

	# Get the B Spline and Farfield B Splines
	BSpline_parameters = get_airfoil_BSpline_parameters()
	Farfield_parameters = get_farfield_BSpline_parameters()

	xiVector = BSpline_parameters["knots"]
	etaVector = CONST_ETA_KNOTS
	P = CONST_P
	Q = CONST_Q

	# Build the control points and weights structure
	num_xi_pts = len(BSpline_parameters["control_points"])
	num_eta_pts = 2

	ControlPoints_and_Weights = []

	for i in range(num_xi_pts):
		col = []
		for j in range(num_eta_pts):
			col.append([None, None, None])
		ControlPoints_and_Weights.append(col)
	
	for i in range(num_xi_pts):
		ControlPoints_and_Weights[i][0] = [BSpline_parameters["control_points"][i][0],
				BSpline_parameters["control_points"][i][1], 1.]
		ControlPoints_and_Weights[i][1] = [Farfield_parameters["control_points"][i][0],
				Farfield_parameters["control_points"][i][1], 1.]

	# Wrap the patch around by closing it. That is, make the last point on the control
	# mesh equal the first one (xi = 1 and xi = -1)
	ControlPoints_and_Weights[-1][0] = ControlPoints_and_Weights[0][0]
	ControlPoints_and_Weights[-1][1] = ControlPoints_and_Weights[0][1]


	# Get the list of optimization control points
	optimization_control_pt_list = get_optimization_pts(ControlPoints_and_Weights)

	patch_parameters = {
		"xiVector" : xiVector,
		"etaVector" : etaVector,
		"ControlPoints_and_Weights" : ControlPoints_and_Weights,
		"P" : P,
		"Q" : Q,
		"Optimization_ControlPoints_and_Weights" : optimization_control_pt_list,
		"area_ref" : 1.0
	}

	return patch_parameters


def test():

	get_patch_information()

	BSpline_parameters = get_airfoil_BSpline_parameters()
	Farfield_parameters = get_farfield_BSpline_parameters()

	
	#plot_airfoil_test()
	plot_spline(BSpline_parameters["Spline_Function"],
				BSpline_parameters["knots"],
				BSpline_parameters["control_points"])

	plot_spline(Farfield_parameters["Spline_Function"],
				Farfield_parameters["knots"],
				Farfield_parameters["control_points"])

	plt.grid()

	#plt.gca().set_xlim(CONST_PlotXRange)
	#plt.gca().set_ylim(CONST_PlotYRange)

	plt.show(block=True)


if __name__ == "__main__":
	test()






