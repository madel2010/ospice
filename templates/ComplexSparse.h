
#ifndef COMPLEX_SPARSE_H
#define COMPLEX_SPARSE_H

#include <iostream>
#include <list>
#include <functional>
#include <algorithm>
#include <vector>
#include <complex>

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>


#include "klu.h"

#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <ctime>

namespace BMatrix{
  


      //This re writes the RHS
      template<>
      inline void Sparse<std::complex< double > >::solve(BMatrix::Dense<std::complex< double > >& RHS, bool do_klu_refactor){
	  
	  this->create_ccs();

	  int Nrhs = RHS.get_number_of_cols();

	  if(this->structure_has_changed && !do_klu_refactor ){	

		if(this->Symbolic){
			klu_free_symbolic (&this->Symbolic, &this->Common);
			this->Symbolic=NULL;
		}

		clock_t start,finish;
		start = clock();

		this->Symbolic = klu_analyze(this->rows, this->Ap, this->Ai, &this->Common) ;
		
		finish = clock();
		this->time_to_do_klu_analyze+= (double(finish)-double(start))/CLOCKS_PER_SEC;
		this->number_of_klu_analyze++;
	  }
	  

	  //Do LU factorization if not done before
	  if(this->structure_has_changed && !do_klu_refactor){
	      
		if(this->Numeric){
			klu_numeric* my_numeric = (klu_numeric*)this->Numeric;
			klu_z_free_numeric (&my_numeric, &this->Common); 
			this->Numeric=NULL;
		}		

	        clock_t start,finish;
		start = clock();

	      	this->Numeric  = klu_z_factor ( this->Ap, this->Ai, reinterpret_cast <double*>(this->Ax), this->Symbolic, &this->Common ) ;

	      	finish = clock();
		this->time_to_do_klu_factor+= (double(finish)-double(start))/CLOCKS_PER_SEC;
		this->number_of_klu_factor++;
		
		//we have done initial LU so set structure_has_changed=false to prevent another full LU if same nz_structure
		this->structure_has_changed = false;

	  }else if(this->values_have_changed || do_klu_refactor) {
		clock_t start,finish;
		start = clock();

		klu_z_refactor ( this->Ap, this->Ai, reinterpret_cast <double*>(this->Ax), this->Symbolic, (klu_numeric*)this->Numeric, &Common ) ;

	      	finish = clock();
		this->time_to_do_klu_refactor+= (double(finish)-double(start))/CLOCKS_PER_SEC;
		this->number_of_klu_refactor++;

		//we have done LU so set values_has_changed=false to prevent another LU if same values
		this->values_have_changed = false;
		this->structure_has_changed = false;
	  }
	  
	  	  
	  //Now do the F/B substitution
	  int result = klu_z_solve(this->Symbolic, (klu_numeric*)this->Numeric, this->rows, Nrhs, (double*)(*RHS), &this->Common);	  
	  if(!result){
		std::cerr<<"Cannot do F/B substitution";
	  }
	  
      }
  

  
};

  

#endif