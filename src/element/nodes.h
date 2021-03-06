/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */

#ifndef DPG__nodes_h__INCLUDED
#define DPG__nodes_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions computing the reference element coordinates (and associated cubature
 *         weights if relevant).
 */

#include <stdbool.h>

/// Container for reference element coordinates and cubature related information.
struct Nodes {
	/** By default, the order of the nodes. For node sets which do not have the correct number of nodes to form a
	 *  polynomial Lagrange basis, this is the order of the cubature strength. */
	int p;

	int node_type; ///< The index associated with the node type.

	struct Matrix_d* rst; ///< The reference element coordinates.

	bool has_weights;   ///< Flag for whether weights are included.
	struct Vector_d* w; ///< The cubature weights.

	bool has_symms;     ///< Flag for whether the symmetries are included.
	struct Vector_i* s; ///< The node symmetries.
};

/// `const` version of \ref Nodes.
struct const_Nodes {
	const int p,         ///< See brief.
	          node_type; ///< See brief.

	const struct const_Matrix_d*const rst; ///< See brief.

	const bool has_weights;              ///< See brief.
	const struct const_Vector_d*const w; ///< See brief.

	const bool has_symms;           ///< See brief.
	const struct const_Vector_i* s; ///< See brief.
};

/** \brief Function pointer to node constructor function.
 *  \param d         The dimension of the nodes.
 *  \param p         Defined in \ref Nodes.
 *  \param node_type Defined in \ref Nodes.
 */
typedef const struct const_Nodes* (*constructor_Nodes_fptr)
	(const int d,
	 const int p,
	 const int node_type
	);

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a \ref Nodes container of tensor-product type.
 *  \return Standard.
 *
 *  The available nodes are:
 *  - GL  : Gauss-Legendre;
 *  - GLL : Gauss-Lobatto-Legendre;
 *  - EQ  : Equally-Spaced (no weights).
 *
 *  The nodes and weights are computed using the Jacobi library extension to the GNU GSL library written by Paulo
 *  Jabardo.
 */
const struct const_Nodes* constructor_const_Nodes_tp
	(const int d,        ///< Defined in \ref constructor_Nodes_fptr.
	 const int p,        ///< Defined in \ref constructor_Nodes_fptr.
	 const int node_type ///< Defined in \ref constructor_Nodes_fptr.
	);

/** \brief Constructor for a \ref Nodes container of simplex type.
 *  \return Standard.
 *
 *  The available nodes are:
 *  - AO  : Alpha-Optimized (\cite Hesthaven2007) (no weights);
 *  - WSH : Williams-Shunn-Ham (\cite Williams2014, \cite Shunn2012);
 *  - WV  : Witherden-Vincent (\cite Witherden2015);
 *  - EQ  : Equally-Spaced (no weights).
 *
 *  The nodes and weights were taken from [pyfr/quadrules][pyfr_web] after conversion to the input format used here.
 *
 *  2d (TRI) WSH nodes have the following [polynomial order, cubature strength] pairs:
 *  - [0,1], [1,2], [2,4], [3,5], [4,7], [5,8], [6,10], [7,12], [8,14].
 *
 *  3d (TET) WSH nodes have the following [polynomial order, cubature strength] pairs:
 *  - [0,1], [1,2], [2,3], [3,5], [4,6], [5,8], [6,9].
 *
 *  <!-- References: -->
 *  [pyfr_web]: http://www.pyfr.org
 */
const struct const_Nodes* constructor_const_Nodes_si
	(const int d,        ///< Defined in \ref constructor_Nodes_fptr.
	 const int p,        ///< Defined in \ref constructor_Nodes_fptr.
	 const int node_type ///< Defined in \ref constructor_Nodes_fptr.
	);

/** \brief Constructor for a \ref Nodes container of pyramid type.
 *  \return Standard.
 *
 *  The available nodes are:
 *  - GL     : Gauss-Legendre         (no weights);
 *  - GLL    : Gauss-Lobatto-Legendre (no weights);
 *  - GLW    : Gauss-Legendre         (with weights);
 *  - GLLW   : Gauss-Lobatto-Legendre (with weights);
 *  - GJW    : Gauss-Jacobi           (with weights);
 *  - WV     : Witherden-Vincent (\cite Witherden2015) (with weights; **integration to lower order than expected**);
 *  - WVHToP : WV HEX nodes transfered to PYR (with weights; **integration to lower order than expected**).
 *
 *  The WV and WVHToP nodes were determined from those of [pyfr/quadrules][pyfr_web]. The GL(W) and GLL(W) nodes are
 *  standard based on the 1d definions and the GJW nodes were generated using alpha = 2.0, beta = 0.0 (cancelling the
 *  (1-c)^2 term resulting from the transformation from PYR to HEX for integration).
 *
 *  The WV nodes integrate exactly (with the expected order) only when there is variation either in the `a` **or** the
 *  `b` directions, but not in both; the nodes cannot exactly integration polynomials on QUAD cross-sections of the
 *  pyramid. Using the WV HEX nodes transferred to the PYR element also does not result in exact integration.
 *
 *  Using the GL/GLL nodes transferred to the PYR element, the cubature strength for exact mass matrix integration is as
 *  expected. After accounting for the added contribution to the `c` term from the weight (w = w_HEX*pow(1-c,2)), GL
 *  nodes of order p+1 are required for exact integration of all terms.
 *
 *  The set of all (symmetrically redundant) nodes is computed from the barycentric coordinates of the node symmetry
 *  groups provided in the input files. Unlike for the simplex elements, barycentric coordinates for the pyramid element
 *  may be negative. This does not result in any limitations, and is simply a result of the last additional condition
 *  when computing the barycentric coordinates of the nodes (chosen for good conditioning).
 *
 *  Given the following conditions:
 *  - [ones]_{Nn x 1}  = [BCoords]_{Nn x Nc} * [ones]_{Nc x 1}           => Partition of unity (1 condition)
 *  - [rst]_{Nn x d}   = [BCoords]_{Nn x Nc} * [rst_(V)ertices]_{Nc x d} => Linear precision   (d conditions)
 *  - [r*s*t]_{Nn x 1} = [BCoords]_{Nn x Nc} * [{r*s*t}_V]_{Nc x 1}      => Arbitrary          (1 condition)
 *
 *  then
 *  	[BCoords] = [ones(Nn,1) rst r*s*t]*inv([ones(Nc,1) rst_V {r*s*t}_V])
 *
 *  where
 *  	cond([ones(Nc,1) rst_V {r*s*t}_V]) = 2.0 (the condition number).
 *
 *  <!-- References: -->
 *  [pyfr_web]: http://www.pyfr.org
 */
const struct const_Nodes* constructor_const_Nodes_pyr
	(const int d,        ///< Defined in \ref constructor_Nodes_fptr.
	 const int p,        ///< Defined in \ref constructor_Nodes_fptr.
	 const int node_type ///< Defined in \ref constructor_Nodes_fptr.
	);

/** \brief Constructor for a \ref Nodes container for the p1 reference element vertices for arbitrary element type.
 *  \return See brief. */
const struct const_Nodes* constructor_const_Nodes_vertices
	(const int d,     ///< Defined in \ref constructor_Nodes_fptr.
	 const int p,     ///< Defined in \ref constructor_Nodes_fptr.
	 const int s_type ///< \ref Element::s_type.
	);

/// \brief Destructor for a \ref Nodes\* container.
void destructor_Nodes
	(struct Nodes* nodes ///< Standard.
	);

/// \brief Destructor for a \ref const_Nodes\* container.
void destructor_const_Nodes
	(const struct const_Nodes*const nodes ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief Get a pointer to the appropriate Nodes constructor function based on the input element super type.
 *  \return See brief. */
constructor_Nodes_fptr get_constructor_Nodes_by_super_type
	(const int s_type ///< \ref Element::s_type.
	);

/// \brief Print a \ref Nodes container to the terminal displaying entries below the tolerance as 0.0.
void print_Nodes_tol
	(const struct Nodes*const nodes, ///< Standard.
	 const double tol                ///< The tolerance.
	);

/// \brief `const` version of \ref print_Nodes_tol.
void print_const_Nodes_tol
	(const struct const_Nodes*const nodes, ///< See brief.
	 const double tol                      ///< See brief.
	);

/// \brief Print a \ref Nodes container to the terminal calling \ref print_Nodes_tol with a default tolerance.
void print_Nodes
	(const struct Nodes*const nodes ///< Standard.
	);

/// \brief `const` version of \ref print_Nodes.
void print_const_Nodes
	(const struct const_Nodes*const nodes ///< See brief.
	);

#endif // DPG__nodes_h__INCLUDED
