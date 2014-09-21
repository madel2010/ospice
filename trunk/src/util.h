//
// C++ Interface: util
//
// Description: 
//
//
// Author: Mina Farhan <madel@assembly.doe.carleton.ca>, (C) 2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "BMatrix.h"
#include <stdlib.h> 

extern "C" int idamax_( const int * n,  const double* x, const int* incx );

//return maximum absolute value of a BMatrix::Dense
double max_abs(const BMatrix::Dense<double> &A){

	int n =A.get_number_of_cols();
	int m = A.get_number_of_rows();
	int size = m*n;
	int idx = 1;

	int index = idamax_(&size , *A, &idx);

	return abs((*A)[index]);
}