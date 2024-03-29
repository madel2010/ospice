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
#ifndef DOUBLE_DENSE_H
#define DOUBLE_DENSE_H

#ifdef __cplusplus

#include <iostream>
#include <stdexcept>
#include <string.h>
#include <MatrixBase.h>


extern "C" void dgetrf_( const int * M,  const int* N, double* A, const int *lda, int* ipiv, int* result );
extern "C" void dgetrs_( const char* TRANS,  const int* N, const int* nrhs, double* A, const int *lda, int* ipiv, double* B, const int* ldb,  int* result, int tlen);
extern "C" double dnrm2_(const int*N , double*A ,const int* incx);
extern "C" void dgemm_( const char* transa,  const char* transb, const int* m, const int* n, const int* k, double* alpha, double* a, const int* lda, double* b, const int* ldb, const double* beta, double* c, const int* ldc, int transalen, int transblen);

namespace BMatrix{
  

/*---------------specialization for double Dense matrices (Dense<double>)-------------*/
 template<> 
 inline Dense<double>::Dense(int m, int n){
	this->rows=m;
	this->cols=n;
	
	//this->data = new double[m*n];
	this->data = new double[m*n];
	bzero(this->data,m*n*sizeof(double));
	
	this->LU_factors = NULL; //we do not need to allocate it unless we are going to invert the matrix
	have_LU_factors = false;
}

template<> 
 inline void Dense<double>::reset(){
	
	bzero(this->data,this->rows*this->cols*sizeof(double));
	
	this->LU_factors = NULL; //we do not need to allocate it unless we are going to invert the matrix
	have_LU_factors = false;
}

template<> 
inline void Dense<double>::create(int m, int n){ 
	this->rows=m;
	this->cols=n;
	
	this->data = new double[m*n];
	bzero(this->data,m*n*sizeof(double));

	this->LU_factors = NULL; //we do not need to allocate it unless we are going to invert the matrix
	have_LU_factors = false;
}
	
template<> 
inline Dense<double>& Dense<double>::operator = (const Dense<double> &A){
		
		//if it is a just created matrix, then create it with the same size of A
		if(!this->data){
		    create(A.rows,A.cols);
		}else if(this->rows!=A.rows || this->cols!=A.cols){
		    throw std::runtime_error("Can not equate two matrices with different sizes");
		}
		
		memcpy(this->data , A.data, this->rows*this->cols*sizeof(double));
	
		LU_factors = NULL; //we are going to allocate it if we are going to invert the matrix
		have_LU_factors = false;
		
		return *this;
}  
	
template<>
inline double Dense<double>::norm(int col){ //col is the column to calculate the norm for (by default col=1) check MatrixBase.h
	      int inc=1;
	      int N = this->rows;
	      
	      double result = dnrm2_(&N, data+(this->rows*col), &inc);
	      return result;
}

//this function finds (this)^{-1}*RHS using LU factorization and F/B subst.
template<> 
inline Dense<double> Dense<double>::solve (double* _RHS , const int Nrhs){
	
	int result;
	int* i_piv = new int[this->rows];
	
	double* RHS = new double[this->cols];
	memcpy(RHS , _RHS , this->cols*sizeof(double));
	
	if(!have_LU_factors){
	    LU_factors = new double[this->rows*this->cols];	
	    memcpy(LU_factors,this->data,this->rows*this->cols*sizeof(double));
	    
	    result = 0;
	    dgetrf_ ( &this->rows, &this->cols, LU_factors, &this->rows, i_piv, &result );	
	    if ( result ){
		    throw std::runtime_error("Cannot find LU factors of a dense matrix");
	    }
	 }

	 result = 0;
	 char tran = 'N';
	 dgetrs_(&tran, &this->rows, &Nrhs, LU_factors, &this->rows, i_piv, RHS, &this->rows,&result , 1);
	 
	 if ( result ){
		  throw std::runtime_error("Cannot Forward/Backward subst for a dense matrix");
	 }
	 
	 delete[] i_piv;
	 
	 return Dense<double>(this->rows, Nrhs , RHS);
}

template<> 
inline Dense<double> Dense<double>::solve (const Dense<double>& _RHS){

	int result;
	int* i_piv = new int[this->rows];

	Dense<double> RHS = _RHS;
	
	if(this->rows==1){
		RHS.data[0] /= this->data[0];
	}else{
		if(!have_LU_factors){
	    		LU_factors = new double[this->rows*this->cols];
	    		memcpy(LU_factors,this->data,this->rows*this->cols*sizeof(double));
	    		result = 0;
	    		dgetrf_ ( &this->rows, &this->cols, LU_factors, &this->rows, i_piv, &result );	
	    		if ( result ){
		    		throw std::runtime_error("Cannot find LU factors of a dense matrix");
	    		}
	 	}

	 	result = 0;
	 	char tran = 'N';
	 	int Nrhs = RHS.cols;
	 	dgetrs_(&tran, &this->rows, &Nrhs, LU_factors, &this->rows, i_piv, RHS.data, &this->rows,&result , 1);
	 
	 	if ( result ){
		  	throw std::runtime_error("Cannot Forward/Backward subst for a dense matrix");
	 	}
	 
	 	delete[] i_piv;
	}

	return RHS;
}


//Note: this function overrides RHS
template<> 
inline DBase<double>* Dense<double>::solve (DBase<double>* _RHS){

	int result;
	int* i_piv = new int[this->rows];

	if(Dense<double>* Casted_RHS=dynamic_cast<Dense<double>* >(_RHS)){
		

		//Dense<double>* RHS = (Casted_RHS);
		
	
		if(this->rows==1){
			Casted_RHS->data[0] /= this->data[0];
		}else{
			if(!have_LU_factors){
	    			LU_factors = new double[this->rows*this->cols];
	    			memcpy(LU_factors,this->data,this->rows*this->cols*sizeof(double));
	    			result = 0;
	    			dgetrf_ ( &this->rows, &this->cols, LU_factors, &this->rows, i_piv, &result );	
	    			if ( result ){
		    			throw std::runtime_error("Cannot find LU factors of a dense matrix");
	    			}
	 		}

	 		result = 0;
	 		char tran = 'N';
	 		int Nrhs = Casted_RHS->cols;
	 		dgetrs_(&tran, &this->rows, &Nrhs, LU_factors, &this->rows, i_piv, Casted_RHS->data, &this->rows,&result , 1);
	 
	 		if ( result ){
		  		throw std::runtime_error("Cannot Forward/Backward subst for a dense matrix");
	 		}
	 
	 		delete[] i_piv;
		}

		return const_cast<Dense<double>*>(Casted_RHS);
	}else{
		throw std::runtime_error("Cannot solve Dense matrix");
	}

	
}

template<> 
inline DBase<double>*  Dense<double>::solve(const DBase<double> &B){
	const Dense<double>* BB = dynamic_cast<const Dense<double>*>(&B);
    
	if(BB){
	    return this->solve( const_cast<Dense<double>*>(BB) );
	}else{
	    throw std::runtime_error("I can only handle solve of Dense+Dense");
	}
}

//// result = inv(A)*B
template<> 
inline Dense<double> Dense<double>::operator / (const Dense<double> &B){
    if(this->cols!=B.rows){
      throw std::runtime_error("Can not divide two matrices with different rows and columns");
    }
    
    //BMatrix::Dense<double> result(B.rows,B.cols);
    
    //result = B;

    return solve(B);
}

template<>
inline Dense<double> Dense<double>::operator * (const Dense<double> &B) const{
	int m = get_number_of_rows();
	int n = B.get_number_of_cols();
	int k = get_number_of_cols();
	int ldb = B.get_number_of_rows();
		
	double Alpha = 1;
	double Beta = 0;
	Dense<double> result(m,n);
		
	//data = alpha*A*B + Beta*data
	dgemm_("N", "N", &m, &n, &k, &Alpha, this->data, &m, *B, &ldb, &Beta, *result, &m,1,1);
	
	return result;
}

/*-----------------------START: BLAS routines--------------------------------*/
//computes this = alpha*A*B + Beta*this
template<>
inline void Dense<double>::multibly_AB_add_this(double Alpha, double Beta, Dense<double>* A, Dense<double>* B){
	
	if(A->get_number_of_cols()!=B->get_number_of_rows() ||
	   A->get_number_of_rows() != this->get_number_of_rows() ||
	   B->get_number_of_cols()!=this->get_number_of_cols()){
		throw std::runtime_error("BlasProduct is trying to multibly two matrices with different cols and rows");
	}
	
	int m = A->get_number_of_rows();
	int n = B->get_number_of_cols();
	int k = A->get_number_of_cols();
	int ldb = B->get_number_of_rows();
	dgemm_("N", "N", &m, &n, &k, &Alpha, **A, &m, **B, &ldb, &Beta, this->data, &(this->rows),1,1);
}	

template<>
inline void Dense<double>::multibly_AB_add_this(double Alpha, double Beta, DBase<double>* A, DBase<double>* B){
	
	if(A->get_number_of_cols()!=B->get_number_of_rows() ||
	   A->get_number_of_rows() != this->get_number_of_rows() ||
	   B->get_number_of_cols()!=this->get_number_of_cols()){
		throw std::runtime_error("BlasProduct is trying to multibly two matrices with different cols and rows");
	}
	
	const Dense<double>* AA = dynamic_cast<const Dense<double>*>(A);
	const Dense<double>* BB = dynamic_cast<const Dense<double>*>(B);

  	if (!BB || !AA) {
    		throw std::runtime_error("I can only handle multibly_AB_add_this of dense");
  	}
  		
	int m = A->get_number_of_rows();
	int n = B->get_number_of_cols();
	int k = A->get_number_of_cols();
	int ldb = B->get_number_of_rows();
	dgemm_("N", "N", &m, &n, &k, &Alpha, **A, &m, **B, &ldb, &Beta, this->data, &(this->rows),1,1);
}
/*-----------------------END: BLAS routines--------------------------------*/

}

#endif
#endif