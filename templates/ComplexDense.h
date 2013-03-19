//
// C++ Interface: MatrixBase
//
// Description: 
//
//
// Author:  <>, (C) 2012
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef COMPLEX_DENSE_H
#define COMPLEX_DENSE_H

#ifdef __cplusplus

#include <iostream>
#include <stdexcept>
#include <string.h>
#include <MatrixBase.h>
#include <complex>


namespace BMatrix{


/*---------------specialization for double Dense matrices (Dense< std::complex<double> >)-------------*/

 template<> 
 inline bool Dense< std::complex<double> >::is_zero(){
	//we only check the diagonals
	for(int i=0; i< this->rows; i++){
		if(real(this->data[i+this->rows*i])!=0 && imag(this->data[i+this->rows*i])!=0){
			return false;
		}
	}
	return true;
}





}

#endif
#endif