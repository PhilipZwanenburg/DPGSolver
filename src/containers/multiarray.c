// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "multiarray.h"

#include <string.h>

#include "macros.h"

#include "matrix.h"
#include "vector.h"

#include "allocators.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for a vector with an associated index.
struct Vector_i_indexed {
	ptrdiff_t index;         ///< The index.
	struct Vector_i* vector; ///< The \ref Vector_i\*.
};

/** Move constructor for a \ref Vector_i_indexed\*\* from a \ref Vector_i\*\*.
 *	\return Standard. */
static struct Vector_i_indexed** constructor_move_Vector_i_indexed
	(const ptrdiff_t size,  ///< The number of elements.
	 struct Vector_i** data ///< The data to be moved.
	);

/// Destructor for a \ref Vector_i_indexed\*\*.
static void destructor_Vector_i_indexed
	(struct Vector_i_indexed** src, ///< Standard.
	 const ptrdiff_t size           ///< The number of elements.
	);

/** \brief Make a local \ref Multiarray_Vector_i\* (dynamic memory).
 *	\return See brief. */
static struct Multiarray_Vector_i* constructor_local_Multiarray_Vector_i_1
	(const int order,            ///< Standard.
	 ptrdiff_t*const extents,    ///< Standard.
	 const bool owns_data,       ///< Standard.
	 struct Vector_i**const data ///< Standard.
	);

/** \brief Make a local \ref Multiarray_d\* (dynamic memory).
 *	\return See brief. */
static struct Multiarray_d* constructor_local_Multiarray_d_1
	(const char layout,       ///< Standard.
	 const int order,         ///< Standard.
	 ptrdiff_t*const extents, ///< Standard.
	 const bool owns_data,    ///< Standard.
	 double*const data        ///< Standard.
	);

/** \brief Comparison function for std::qsort between `struct Vector_i_indexed**` `a` and `b`.
 *	\return The lexicographical comparison of `a` and `b`.
 *
 *	\note Input Vectors must be have sorted data.
 */
static int cmp_Vector_i_indexed
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/** \brief Reorder a \ref Multiarray_Vector_i based on the provided ordering.
 *	\warning This is not currently done in place.
 */
static void reorder_Multiarray_Vector_i
	(struct Multiarray_Vector_i*const a, ///< Standard.
	 const int*const ordering    ///< The ordering.
	);

/// \brief Compute the total number of entries in a \ref Multiarray_Vector_i\*.
static ptrdiff_t compute_total_entries
	(const struct Multiarray_Vector_i*const src ///< Standard.
	);

/// \brief Print the extents of the Multiarray.
void print_Multiarray_extents
	(const int order,              ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents ///< Defined in \ref Multiarray_d.
	);

/// \brief Print the counter for the indices of order > 2 when printing sub-Matrices of the Multiarray.
void print_Multiarray_counter
	(const int order,              ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const counter ///< The counter for the indices of order > 2.
	);

/**	\brief Compute the index of the data of the sub-Matrix in a Multiarray based on the counter.
 *	\return See brief. */
static ptrdiff_t compute_ind_sub_matrix
	(const int order,               ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents, ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const counter  ///< The counter for the indices of order > 2.
	);

///	\brief Increment the counter for the indicies of order > 2 by 1.
static void increment_counter
	(const int order,               ///< Defined in \ref Multiarray_d.
	 const ptrdiff_t*const extents, ///< Defined in \ref Multiarray_d.
	 ptrdiff_t*const counter        ///< The counter for the indices of order > 2.
	);

// Constructor/Destructor functions ********************************************************************************* //

struct Multiarray_d* constructor_default_Multiarray_d ()
{
	const int order = 1;
	ptrdiff_t* extents = calloc(order , sizeof *extents); // keep

	return constructor_local_Multiarray_d_1('C',order,extents,true,NULL);
}

struct Multiarray_d* constructor_move_Multiarray_d_d (const char layout, double*const data, const int order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	ptrdiff_t* extents = allocate_and_set_extents(order,ap); // keep
	va_end(ap);

	return constructor_local_Multiarray_d_1(layout,order,extents,false,data);
}

struct Multiarray_Vector_i* constructor_empty_Multiarray_Vector_i (const bool alloc_V, const int order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	ptrdiff_t* extents = allocate_and_set_extents(order,ap); // keep
	va_end(ap);

	struct Vector_i** data = NULL;
	if (alloc_V)
		data = constructor_default_Vector_i_2(compute_size(order,extents)); // keep
	else
		data = malloc(compute_size(order,extents) * sizeof *data); // keep

	return constructor_local_Multiarray_Vector_i_1(order,extents,true,data);
}

struct Multiarray_Vector_i* constructor_copy_Multiarray_Vector_i_i
	(const int* data_V, const int*const ext_V, const int order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	ptrdiff_t* extents = allocate_and_set_extents(order,ap); // free (reallocated in constructor_empty)
	va_end(ap);

	ptrdiff_t extents_l[order];
	for (ptrdiff_t i = 0; i < order; ++i)
		extents_l[i] = extents[i];
	free(extents);

	struct Multiarray_Vector_i* dest;
	switch (order) {
	case 1:
		dest = constructor_empty_Multiarray_Vector_i(true,order,extents_l[0]);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	set_Multiarray_Vector_i_i(dest,data_V,ext_V);

	return dest;
}

void const_constructor_move_Multiarray_d (const struct const_Multiarray_d*const* dest, struct Multiarray_d* src)
{
	*(struct const_Multiarray_d**) dest = (struct const_Multiarray_d*) src;
}

void const_constructor_move_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const* dest, struct Multiarray_Vector_i* src)
{
	*(struct const_Multiarray_Vector_i**) dest = (struct const_Multiarray_Vector_i*) src;
}

void destructor_Multiarray_d (struct Multiarray_d* a)
{
	free(a->extents);
	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_Multiarray_Vector_i (struct Multiarray_Vector_i* a)
{
	destructor_Vector_i_2(a->data,compute_size(a->order,a->extents),a->owns_data);
	free(a->extents);
	free(a);
}

// Helper functions ************************************************************************************************* //

ptrdiff_t* allocate_and_set_extents (const int order, va_list ap)
{
	ptrdiff_t* extents = mallocator(PTRDIFF_T,1,order); // returned

	for (ptrdiff_t i = 0; i < order; i++)
		extents[i] = va_arg(ap,ptrdiff_t);

	return extents;
}

ptrdiff_t compute_size (const int order, const ptrdiff_t *const extents)
{
	ptrdiff_t size = 1;
	for (ptrdiff_t i = 0; i < order; i++)
		size *= extents[i];

	return size;
}

void set_Multiarray_Vector_i_i
	(struct Multiarray_Vector_i* a, const int* data_V, const int*const ext_V)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; i++) {
		resize_Vector_i(a->data[i],ext_V[i]);
		for (ptrdiff_t j = 0; j < ext_V[i]; j++)
			a->data[i]->data[j] = *data_V++;
	}
}

struct Vector_i* sort_Multiarray_Vector_i (struct Multiarray_Vector_i* a, const bool return_indices)
{
	// sort the individual Vector entries
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		sort_Vector_i(a->data[i]);

	// sort the Vectors
	int* ordering_i = NULL;
	if (!return_indices) {
		qsort(a->data,size,sizeof(a->data[0]),cmp_Vector_i);
	} else {
		struct Vector_i_indexed** a_indexed = constructor_move_Vector_i_indexed(size,a->data); // destructed

		qsort(a_indexed,size,sizeof(a_indexed[0]),cmp_Vector_i_indexed);

		ordering_i = malloc(size * sizeof *ordering_i); // keep
		for (ptrdiff_t i = 0; i < size; ++i)
			ordering_i[i] = a_indexed[i]->index;

		reorder_Multiarray_Vector_i(a,ordering_i);

		destructor_Vector_i_indexed(a_indexed,size);
	}

	struct Vector_i* ordering = constructor_move_Vector_i_i(size,true,ordering_i); // returned

	return ordering;
}

struct Vector_i* collapse_Multiarray_Vector_i (const struct Multiarray_Vector_i*const src)
{
	const ptrdiff_t n_entries = compute_total_entries(src);

	int*const data = malloc(n_entries * sizeof *data); // keep

	ptrdiff_t ind_d = 0;
	const ptrdiff_t size = compute_size(src->order,src->extents);
	for (ptrdiff_t i = 0; i < size; ++i) {
		struct Vector_i*const src_curr = src->data[i];
		const ptrdiff_t size_V = src_curr->ext_0;
		for (ptrdiff_t j = 0; j < size_V; ++j) {
			data[ind_d] = src_curr->data[j];
			++ind_d;
		}
	}

	struct Vector_i*const dest = calloc(1,sizeof *dest); // returned
	dest->ext_0     = n_entries;
	dest->owns_data = true;
	dest->data      = data;

	return dest;
}

// Printing functions *********************************************************************************************** //

void print_Multiarray_Vector_i (const struct Multiarray_Vector_i*const a)
{
	print_Multiarray_extents(a->order,a->extents);

	const ptrdiff_t size = compute_size(a->order,a->extents);

	const int order = a->order;
	switch (order) {
	case 1:
		for (ptrdiff_t i = 0; i < size; i++)
			print_Vector_i(a->data[i]);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	printf("\n");
}

void print_const_Multiarray_Vector_i (const struct const_Multiarray_Vector_i*const a)
{
	struct Multiarray_Vector_i* local =
		constructor_local_Multiarray_Vector_i_1(a->order,(ptrdiff_t*)a->extents,false,(struct Vector_i**)a->data);
	print_Multiarray_Vector_i(local);
	free(local);
}

void print_Multiarray_d (const struct Multiarray_d*const a, const double tol)
{
	const int order               = a->order;
	const ptrdiff_t*const extents = a->extents;

	print_Multiarray_extents(order,extents);

	if (order == 1) {
		struct Vector_d* a_V = constructor_move_Vector_d_d(extents[0],false,a->data); // destructed
		print_Vector_d(a_V,tol);
		destructor_Vector_d(a_V);
	} else if (order == 2) {
		struct Matrix_d* a_M = constructor_move_Matrix_d_d(a->layout,extents[0],extents[1],false,a->data); // destructed
		print_Matrix_d(a_M,tol);
		destructor_Matrix_d(a_M);
	} else if (order > 2) {
		struct Matrix_d* a_M = constructor_move_Matrix_d_d(a->layout,extents[0],extents[1],false,NULL); // destructed;

		const ptrdiff_t size_tail = compute_size(order-2,&extents[2]);

		ptrdiff_t counter[size_tail];
		for (ptrdiff_t i = 0; i < size_tail; ++i)
			counter[i] = 0;

		for (ptrdiff_t i = 0; i < size_tail; ++i) {
			print_Multiarray_counter(order,counter);

			const ptrdiff_t ind_M = compute_ind_sub_matrix(order,extents,counter);
			a_M->data = &a->data[ind_M];
			print_Matrix_d(a_M,tol);

			increment_counter(order,extents,counter);
		}
		destructor_Matrix_d(a_M);
	} else {
		EXIT_UNSUPPORTED;
	}
	printf("\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Vector_i_indexed** constructor_move_Vector_i_indexed
	(const ptrdiff_t size, struct Vector_i** data)
{
	struct Vector_i_indexed** dest = malloc(size * sizeof *dest); // returned

	for (ptrdiff_t i = 0; i < size; i++) {
		dest[i] = calloc(1,sizeof *dest[i]); // keep
		dest[i]->index  = i;
		dest[i]->vector = data[i];
	}
	return dest;
}

static void destructor_Vector_i_indexed (struct Vector_i_indexed** src, const ptrdiff_t size)
{
	for (ptrdiff_t i = 0; i < size; ++i)
		free(src[i]);
	free(src);
}

static struct Multiarray_d* constructor_local_Multiarray_d_1
	(const char layout, const int order, ptrdiff_t*const extents, const bool owns_data, double*const data)
{
	struct Multiarray_d* dest = calloc(1,sizeof *dest); // returned

	dest->layout    = layout;
	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

static struct Multiarray_Vector_i* constructor_local_Multiarray_Vector_i_1
	(const int order, ptrdiff_t*const extents, const bool owns_data, struct Vector_i**const data)
{
	struct Multiarray_Vector_i* dest = calloc(1,sizeof *dest); // returned

	dest->order     = order;
	dest->extents   = extents;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

static void reorder_Multiarray_Vector_i (struct Multiarray_Vector_i*const a, const int*const ordering)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	struct Vector_i* b[size];
	for (ptrdiff_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (ptrdiff_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

static int cmp_Vector_i_indexed (const void *a, const void *b)
{
	const struct Vector_i_indexed*const*const ia = (const struct Vector_i_indexed*const*const) a,
	                             *const*const ib = (const struct Vector_i_indexed*const*const) b;

	const ptrdiff_t size_a = (*ia)->vector->ext_0,
	                size_b = (*ib)->vector->ext_0;

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const int*const data_a = (*ia)->vector->data,
	         *const data_b = (*ib)->vector->data;

	for (ptrdiff_t i = 0; i < size_a; ++i) {
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
	}
	return 0;
}

static ptrdiff_t compute_total_entries (const struct Multiarray_Vector_i*const src)
{
	const ptrdiff_t size = compute_size(src->order,src->extents);

	ptrdiff_t count = 0;
	for (ptrdiff_t i = 0; i < size; ++i)
		count += src->data[i]->ext_0;

	return count;
}

void print_Multiarray_extents (const int order, const ptrdiff_t*const extents)
{
	printf("Multi-array extents: {");
	for (ptrdiff_t i = 0; i < order; i++)
		printf(" %td,",extents[i]);
	printf(" }\n\n");
}

void print_Multiarray_counter (const int order, const ptrdiff_t*const counter)
{
	printf("{:,:");
	for (int i = 0; i < order-2; ++i)
		printf(",%td",counter[i]);
	printf("}\n");
}

static ptrdiff_t compute_ind_sub_matrix (const int order, const ptrdiff_t*const extents, const ptrdiff_t*const counter)
{
	const ptrdiff_t*const extents_tail = &extents[2];

	ptrdiff_t base = extents[0]*extents[1];
	ptrdiff_t ind_sub_m = 0;
	for (int i = 0; i < order-2; ++i) {
		ind_sub_m += base*counter[i];
		base *= extents_tail[i];
	}
	return ind_sub_m;
}

static void increment_counter (const int order, const ptrdiff_t*const extents, ptrdiff_t*const counter)
{
	const ptrdiff_t*const extents_tail = &extents[2];

	for (int i = 0; i < order-2; ++i) {
		if (counter[i] == extents_tail[i]) {
			counter[i] = 0;
			continue;
		}
		++counter[i];
		return;
	}
	EXIT_ERROR("Counter greater than extents.");
}
