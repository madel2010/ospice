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
#ifndef MATRIX_WRAP_H
#define MATRIX_WRAP_H

#ifdef __cplusplus

#include <iostream>
#include <stdexcept>
#include <string.h>
#include <MatrixBase.h>


namespace BMatrix{


template <class T>
class MWrap{

private:
	DBase<T> *A;

public:
	double scalar; //this will be needed in the REAL() operator in KLU lib. It is needed for the scaling matrix RS........check klu_factor.cc around line 440
		       
	static bool start_stealing_data; //this is used in KLU_SOLVE to steal the matrix in the case of operator= in order to avoid un-needed copying the data.

	MWrap(DBase<T>*_A) : A(_A) { }
	MWrap() : A(NULL) { }
  	~MWrap() { }
	
	void clear(){
	      //if(A){
	      //	  delete A;
	      //	  A = NULL;
	      //}
	}
	
	void exchange_A(MWrap<T>& B){
	    DBase<T>* _A = this->A;
	    A = B.A;
	    B.A = _A;
	}
	
	//derefernce operator (*A)
	virtual DBase<T>* operator* () const {return A;}

	virtual MWrap<T>& operator +=(const MWrap<T> &B){
		(*A)+=(**B);
		return *this;
	}

	virtual MWrap<T>& operator -=(const MWrap<T>  &B){
		(*A)-=(**B);
		return *this;
	}

	virtual MWrap<T>& operator *=(const MWrap<T>  &B){
		(*A)*=(**B);
		return *this;
	}

	//virtual MWrap<T>& operator /=(MWrap<T>  &B){
		//(*A)/=(**B);
		//return *this;
	//}

	/*virtual MWrap<T>& operator /=(T val){

		Dense<T>* BB = dynamic_cast<Dense<T>*>(A);
		if(BB){
			(*BB)/=val;
		}
		return *this;
	}*/

	virtual MWrap<T>& operator =(T val){
		
		//first check that Mwrap already has a matrix. If not, then do nothing
		if(A){
		  (*A)=val;
		}
		return *this;
	};

	virtual MWrap<T>& operator =(const MWrap<T>  &B){

		if(! MWrap<T>::start_stealing_data){
			if(!A){
		    		if (dynamic_cast<Dense<T>*>(B.A)){ 
					A = new Dense<T>(B.A->get_number_of_rows(),B.A->get_number_of_cols());
		    		}
			}
			(*A) = (**B);
		}else{
			A = B.A;
		}
		return *this;
	};
	
	virtual MWrap<T>& operator=(DBase<T>* B){
		//if(!A){
		//    if (dynamic_cast<Dense<T>*>(B)){ 
		//	A = new Dense<T>;
		//    }
		    
		//}
		
		A = B;
		return *this;
		
		
	};

	virtual MWrap<T> operator *(const MWrap<T>  &B) const{
		return MWrap<T>( *A * *(*B) );
	};

	virtual MWrap<T> operator -(const MWrap<T>  &B) const{
		return MWrap<T>( *A - *(*B) );
	};

	virtual MWrap<T> operator +(const MWrap<T>  &B) const{
		return MWrap<T>( *A + *(*B) );
	};
	
	virtual MWrap<T> operator /(T val) const{
		return MWrap<T>( A->scale(&val) );
	};
	
	virtual MWrap<T>& operator /=(T val){
		*A /= val ;
		return *this;
	};

	virtual MWrap<T> solve(const MWrap<T>  &B) const{
		return MWrap<T>( A->solve(*(*B)) );
	};

	virtual double det(){
		//return A->det();
		double det = 1;
		int rows = A->get_number_of_rows();
		for(int i=0; i< rows; i++){
		    det *= A->get(i,i);
		}
	      
		return det;
	}
	
	virtual double is_zero(){
		return A->is_zero();
	}

};

template <class T> bool MWrap<T>::start_stealing_data = false;

};

#endif
#endif