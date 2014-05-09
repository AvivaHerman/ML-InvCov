/*---------------------------------------------------------------------------
Copyright (2014): Eran Treister, Aviva Herman and Irad Yavneh. 
This code is distributed under the terms of the GNU General Public
License 2.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose."
---------------------------------------------------------------------------*/
#include "sort.h"

static void q_sort_ascend(size_t* indices, size_t left, size_t right);
static void q_sort_descend(double* numbers, size_t* indices, size_t left, size_t right);
static void q_sort_ascend_numbers(double* numbers, size_t* indices, size_t left, size_t right);
static void negate_numbers(double* numbers, size_t array_size);

void quickSort_ascend(double* numbers, size_t* indices, size_t array_size)
{
	int j;
	if (!numbers) {
		// mexPrintf("Im here2\n");
		// if (indices!=0){
		q_sort_ascend(indices, 0, array_size - 1);
		// }
		
		// for (j=0; j<array_size; ++j) {
			// mexPrintf("%d, ", indices[j]);
		// }
	} else {
		// mexPrintf("Im here\n");
		// q_sort_ascend(indices, 0, array_size - 1);
		q_sort_ascend_numbers(numbers, indices, 0, array_size - 1);
		// negate_numbers(numbers, array_size);
		// q_sort_descend(numbers, indices, 0, array_size - 1);
		// negate_numbers(numbers, array_size);
	}
}


void quickSort_descend(double* numbers, size_t* indices, size_t array_size)
{
    q_sort_descend(numbers, indices, 0, array_size - 1);
}


static void q_sort_ascend(size_t* indices, size_t left, size_t right)
{
	size_t pivot_idx, l_hold, r_hold;
  
	l_hold = left;
	r_hold = right;
	pivot_idx = indices[left];

	while (left < right) {
		while ((indices[right] >= pivot_idx) && (left < right))
		  right--;
		  
		if (left != right) {
		  indices[left] = indices[right];
		  left++;
		}

		while ((indices[left] <= pivot_idx) && (left < right))
		  left++;

		if (left != right) {
		  indices[right] = indices[left];
		  right--;
		}
	}

	indices[left] = pivot_idx;
	pivot_idx = left;
	left = l_hold;
	right = r_hold;

	if (left < pivot_idx)
		q_sort_ascend(indices, left, pivot_idx - 1);

	if (right > pivot_idx)
		q_sort_ascend(indices, pivot_idx + 1, right);
}


static void q_sort_ascend_numbers(double* numbers, size_t* indices, size_t left, size_t right)
{
	size_t pivot_idx, l_hold, r_hold;
	double pivot;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];
	pivot_idx = indices[left];

	while (left < right) {
		while ((numbers[right] >= pivot) && (left < right))
		  right--;

		if (left != right) {
		  numbers[left] = numbers[right];
		  indices[left] = indices[right];
		  left++;
		}

		while ((numbers[left] <= pivot) && (left < right))
			left++;

		if (left != right) {
		  numbers[right] = numbers[left];
		  indices[right] = indices[left];
		  right--;
		}
	}

	numbers[left] = pivot;
	indices[left] = pivot_idx;
	pivot_idx = left;
	left = l_hold;
	right = r_hold;

	if (left < pivot_idx)
		q_sort_ascend_numbers(numbers, indices, left, pivot_idx - 1);

	if (right > pivot_idx)
		q_sort_ascend_numbers(numbers, indices, pivot_idx + 1, right);
}


static void q_sort_descend(double* numbers, size_t* indices, size_t left, size_t right)
{
	size_t pivot_idx, l_hold, r_hold;
	double pivot;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];
	pivot_idx = indices[left];

	while (left < right) {
		while ((numbers[right] <= pivot) && (left < right))
		  right--;

		if (left != right) {
		  numbers[left] = numbers[right];
		  indices[left] = indices[right];
		  left++;
		}

		while ((numbers[left] >= pivot) && (left < right))
			left++;

		if (left != right) {
		  numbers[right] = numbers[left];
		  indices[right] = indices[left];
		  right--;
		}
	}

	numbers[left] = pivot;
	indices[left] = pivot_idx;
	pivot_idx = left;
	left = l_hold;
	right = r_hold;

	if (left < pivot_idx)
		q_sort_descend(numbers, indices, left, pivot_idx - 1);

	if (right > pivot_idx)
		q_sort_descend(numbers, indices, pivot_idx + 1, right);
}


static void negate_numbers(double* numbers, size_t array_size){
	int i;
	for (i = 0; i < array_size; ++i) {
		numbers[i] = -numbers[i];
	}
}
