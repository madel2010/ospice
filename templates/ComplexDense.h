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
inline Dense<std::complex<double> >& Dense<std::complex<double> >::operator = (const Dense<std::complex<double> > &A){
		
		//if it is a just created matrix, then create it with the same size of A
		if(!this->data){
		    create(A.rows,A.cols);
		}else if(this->rows!=A.rows || this->cols!=A.cols){
		    throw std::runtime_error("Can not equate two matrices with different sizes");
		}
		
		memcpy(this->data , A.data, this->rows*this->cols*sizeof(std::complex<double>));
	
		LU_factors = NULL; //we are going to allocate it if we are going to invert the matrix
		have_LU_factors = false;
		
		return *this;
}  

 template<> 
 inline Dense<std::complex<double> >::Dense(int m, int n){
	this->rows=m;
	this->cols=n;
	
	//this->data = new double[m*n];
	this->data = new std::complex<double>[m*n];
	bzero(this->data,m*n*sizeof(std::complex<double>));
	
	this->LU_factors = NULL; //we do not need to allocate it unless we are going to invert the matrix
	have_LU_factors = false;
}



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