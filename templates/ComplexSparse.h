
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
  

template<>
class Sparse<std::complex<double> >: public SBase<std::complex<double> >{

private:
      int nnz; //this will save the number of nonzeros
      bool ccs_created;
      
      //For the CCS 
      int* Ap;  //the pointer to columns position
      int* Ai; //the rows data
      std::complex<double>* Ax; 		//The values
      
      //For KLU
      klu_symbolic* Symbolic;
      klu_common Common;
      klu_numeric* Numeric;
      bool structure_has_changed;
      bool values_have_changed;
      

      //the following numbers are for storing the time for doing sparse order and LU factorization.
      double time_to_do_klu_analyze;
      double time_to_do_klu_factor;
      double time_to_do_klu_refactor;
      int number_of_klu_analyze;
      int number_of_klu_factor;
      int number_of_klu_refactor;

      struct SparseElement{
	  int row;
	  std::complex<double> value;
      };
      
      struct FindRow: public std::binary_function< SparseElement, int, bool > {
		    bool operator () ( const SparseElement & element, const int &row ) const {
		    return element.row == row;
		}
      };
      std::list<SparseElement> *cols_lists;
      
            
      //Create the CCS formal from the linked list
      void create_ccs(){
	  if(ccs_created==true){
		return;  //we do not need to redo it if we already have found the css
	  }  
	  
	  Ap[0] = 0; //Ap[0] is always zero
	  
	  Ai = new  int[nnz];
	  Ax = new  std::complex<double>[nnz]; 
	  
	  
	  //typedef typename std::list<SparseElement>::iterator row_iter;
	  std::list<SparseElement>::iterator row_iter;
	  int k = 0;
	  

	  for(int i=0; i< this->cols; i++){
	      Ap[i+1] = Ap[i] + this->cols_lists[i].size();
	      
	      
	      for(row_iter = this->cols_lists[i].begin(); row_iter!=this->cols_lists[i].end(); row_iter++){
		  Ai[k] = row_iter->row;
		  Ax[k] = row_iter->value; 
		  
		  k++;
	      }
	      
	      
	  }  
	  
	  ccs_created = true;
	  
      }
      
      
      void create_MWrap(){}
      
public:
      Sparse(){
	  this->rows=0; //Number of rows
	  this->cols=0; //Number of cols
	  
	  Ap = NULL;
	  Ai = NULL;
	  Ax = NULL;
	  Symbolic = NULL;
	  Numeric = NULL;

	  klu_defaults(&Common);
	  Common.scale=0;

	  ccs_created = false; 
	  structure_has_changed = true;
	  values_have_changed = true;

	  nnz = 0;

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
      }
      
      Sparse(int m, int n){
	  this->rows=m; //Number of rows
	  this->cols=n; //Number of cols
	  
	  cols_lists = new std::list<SparseElement>[n];
	  
	  Ap = new int[n+1]; //We know that Ap is always (columns_no+1)
	  Ai = NULL; //we still do not know the rows
	  Ax = NULL; //we still do not know the values
	  Symbolic = NULL;
	  Numeric = NULL;

	  klu_defaults(&Common);
	  Common.scale=0;

	  ccs_created = false;
	  structure_has_changed = true;
	  values_have_changed = true;

	  nnz = 0;

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
      }
      
      //copy constructor
      Sparse(const Sparse<std::complex<double> >&A){
	  create(A.rows, A.cols);
	  
	  nnz = A.nnz;
	  this->rows=A.rows; //Number of rows
	  this->cols=A.cols; //Number of cols
	  
	  
	  
	  if(A.cols > 0){
	      for(int i=0 ; i< A.cols; i++){
		  cols_lists[i] = A.cols_lists[i];
	      }
	  }  
	  
	  //TODO copy the CCS structure as well
	  ccs_created = false;
	
	  structure_has_changed = true;
	  values_have_changed = true;

	  Symbolic = NULL;
	  Numeric = NULL;

	  klu_defaults(&Common);
	  Common.scale=0;

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
      }

       //convert a double matrix to complex matrix
       Sparse(const Sparse<double> &A){
	  create(A.rows, A.cols);
	  
	  nnz = A.nnz;
	  this->rows=A.rows; //Number of rows
	  this->cols=A.cols; //Number of cols
	  
	  
	  

	  if(A.cols > 0){
	      std::list<Sparse<double>::SparseElement>::iterator row_iter;
	      for(int i=0 ; i< A.cols; i++){
		        for(row_iter = (A.cols_lists[i]).begin(); row_iter!=(A.cols_lists[i]).end(); row_iter++){
		  		SparseElement new_element; //create a new element structure
		  		new_element.row = row_iter->row;  //add the row to the new element structure
		  		new_element.value = row_iter->value; ; //add the value to the new element structure
			 
		  		cols_lists[i].push_back(new_element); 	
	      		}
	      }
	  }  
	  
	  //TODO copy the CCS structure as well
	  ccs_created = false;
	
	  structure_has_changed = true;
	  values_have_changed = true;

	  Symbolic = NULL;
	  Numeric = NULL;
	  klu_defaults(&Common);
	  Common.scale=0;

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
      }
      
	Sparse<std::complex<double> >* clone(){
		return new Sparse<std::complex<double> >(*this);
	}
	
	void create(int m, int n){ 
		this->rows=m; //Number of rows
		this->cols=n; //Number of cols
	  
		cols_lists = new std::list<SparseElement>[n];
	  
		Ap = new int[n+1]; //We know that Ap is always (columns_no+1)
		Ai = NULL; //we still do not know the rows
		Ax = NULL; //we still do not know the values
		Symbolic = NULL;
        	Numeric = NULL;

		klu_defaults(&Common);
	  	Common.scale=0;

		ccs_created = false;
		structure_has_changed = true;
	  	values_have_changed = true;
	  	
		nnz = 0;

		//the time report data
	  	time_to_do_klu_analyze = 0;
      	  	time_to_do_klu_factor = 0;
	  	time_to_do_klu_refactor = 0;
          	number_of_klu_analyze = 0;
          	number_of_klu_factor = 0;
	  	number_of_klu_refactor = 0;
	}
	
      ~Sparse(){
	if(Ap) delete[] Ap;
	Ap=NULL;
	if(Ai) delete[] Ai;
	Ai = NULL;
	if(Ax) delete[] Ax;	

	klu_free_symbolic (&Symbolic, &Common);
	Symbolic=NULL;
        klu_z_free_numeric (&Numeric, &Common); 
	Numeric=NULL;

      }

      Sparse<std::complex<double> >& operator=(std::complex<double> val){}
      
      SBase<std::complex<double> >* scale(std::complex<double>* val)const {}
	
      
      
      Sparse<std::complex<double> >& operator=(const Sparse<std::complex<double> >&A){
	  nnz = A.nnz;
	  this->rows=A.rows; //Number of rows
	  this->cols=A.cols; //Number of cols
	  
	  if(A.cols > 0){
	      for(int i=0 ; i< A.cols; i++){
		  cols_lists[i] = A.cols_lists[i];
	      }
	  }  
	  
	  //TODO copy the CCS structure as well
	  ccs_created = false; //for now we just make the new object create the CCS
	  structure_has_changed = true;
	  values_have_changed = true;

	 return *this;
      }
      
      //Get a value at row m and column n
      std::complex<double> get(int m, int n) const{
	  std::complex<double> result;
	  
	  //search column for the row
	  //typename std::list<SparseElement>::iterator row_iterator;
	  std::list<SparseElement>::iterator row_iterator;

	 row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	 
	 if(row_iterator!=cols_lists[n].end()){ //we have found the row
	    result = row_iterator->value;
	    
	 }else{
	    result = 0;
	 }
	 
	 return result;
	  
      }
      
      std::complex<double> det(){

      }

      bool is_zero(){

      }

      //add a  value in row m and column n to the existing value
      void add_to_entry(int m, int n, std::complex<double> value){
	  
	    if(real(value)==0.0 && imag(value)==0.0) return;

	   //search column for the row in case we have added the row already
	  //typename std::list<SparseElement>::iterator row_iterator;
	  std::list<SparseElement>::iterator row_iterator;

	 row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	 
	 
	 if(row_iterator!=cols_lists[n].end()){ //we have found a value already
	    row_iterator->value += value;
	    values_have_changed = true;

	 }else{  //we have no value at this row. We have to add the row first
	      row_iterator = cols_lists[n].begin();
	      
	      //if no row has been already added, then just add the row
	      if(row_iterator==cols_lists[n].end()){
		  SparseElement new_element; //create a new element structure
		  new_element.row = m;  //add the row to the new element structure
		  new_element.value = value; //add the value to the new element structure
			 
		  cols_lists[n].insert(row_iterator , new_element); 	 
		  nnz++; //increase the number of non zeros in the matrix
		  
	      }else{
		  bool found_larger_row = false;
		  while(row_iterator!=cols_lists[n].end()){
		      if(row_iterator->row > m){  //we have found the first row greater than the required row
		    
			    SparseElement new_element; //create a new element structure
			    new_element.row = m;  //add the row to the new element structure
			    new_element.value = value; //add the value to the new element structure
			 
			    cols_lists[n].insert(row_iterator , new_element); 
			 
			    nnz++; //increase the number of non zeros in the matrix
			 
			    found_larger_row = true;
			    break;
		      }
		      row_iterator++;
		  }

		  //if we didnot find any larger row, then it means that all the rows are smaller, we have to add the new entry at the end
		  if(!found_larger_row){
			  SparseElement new_element; //create a new element structure
			  new_element.row = m;  //add the row to the new element structure
			  new_element.value = value; //add the value to the new element structure
			 
			  cols_lists[n].push_back(new_element); 
			 
			  nnz++; //increase the number of non zeros in the matrix
		  }
	      }

	      
	      structure_has_changed = true;
	 }
	 
	 if(ccs_created){
		if(Ai) delete[] Ai;
		Ai = NULL;
		if(Ax) delete[] Ax;
		Ax = NULL;
	 }
	 ccs_created = false;
	 
      }
      
      //put value in row m and column n
      void put(int m, int n, std::complex<double> value){
	  
	  if(real(value)==0.0 && imag(value)==0.0) return;
	  
	  //search column for the row in case we have added the row already
	  //typename std::list<SparseElement>::iterator row_iterator;
	  std::list<SparseElement>::iterator row_iterator;

	  row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	 
	 
	  if(row_iterator!=cols_lists[n].end()){ //we have found a value already
	     row_iterator->value = value;
	     values_have_changed = true;
	  }else{  //we have no value at this row. We have to add the row first
	       row_iterator = cols_lists[n].begin();
	      
	       //if no row has been already added, then just add the row
	       if(row_iterator==cols_lists[n].end()){
		   SparseElement new_element; //create a new element structure
		   new_element.row = m;  //add the row to the new element structure
		   new_element.value = value; //add the value to the new element structure
			 
		   cols_lists[n].insert(row_iterator , new_element); 	 
		   nnz++; //increase the number of non zeros in the matrix
		  
	       }else{
		   bool found_larger_row = false;
		   while(row_iterator!=cols_lists[n].end()){
		       if(row_iterator->row > m){  //we have found the first row greater than the required row
		    
			     SparseElement new_element; //create a new element structure
			     new_element.row = m;  //add the row to the new element structure
			     new_element.value = value; //add the value to the new element structure
			 
			     cols_lists[n].insert(row_iterator , new_element); 
			 
			     nnz++; //increase the number of non zeros in the matrix
			 
			     found_larger_row = true;
			     break;
		       }
		       row_iterator++;
		   }

		   //if we didnot find any larger row, then it means that all the rows are smaller, we have to add the new entry at the end
		   if(!found_larger_row){
			  SparseElement new_element; //create a new element structure
			  new_element.row = m;  //add the row to the new element structure
			  new_element.value = value; //add the value to the new element structure
			 
			  cols_lists[n].push_back(new_element); 
			 
			  nnz++; //increase the number of non zeros in the matrix
		  }
	       }
	
	       structure_has_changed = true;
	  }
	  
	  if(ccs_created){
		if(Ai) delete[] Ai;
		Ai = NULL;
		if(Ax) delete[] Ax;
		Ax = NULL;
	  }
	  ccs_created = false; 
      }


 
      //get the number of non zeros entries
      int get_nnz(){ return nnz; }
      
      //This re writes the RHS
      void solve(BMatrix::Dense<std::complex<double> >& RHS, bool do_klu_refactor=false){
	  
	  create_ccs();

	  int Nrhs = RHS.get_number_of_cols();

	  if(structure_has_changed && !do_klu_refactor ){	
		clock_t start,finish;
		start = clock();

		Symbolic = klu_analyze(this->rows, Ap, Ai, &Common) ;
		
		finish = clock();
		time_to_do_klu_analyze+= (double(finish)-double(start))/CLOCKS_PER_SEC;
		number_of_klu_analyze++;
	  }
	  

	  //Do LU factorization if not done before
	  if(structure_has_changed && !do_klu_refactor){
	      
	        clock_t start,finish;
		start = clock();

	      	Numeric  = klu_z_factor ( Ap, Ai, reinterpret_cast <double*>(Ax), Symbolic, &Common ) ;

	      	finish = clock();
		time_to_do_klu_factor+= (double(finish)-double(start))/CLOCKS_PER_SEC;
		number_of_klu_factor++;
		
		//we have done initial LU so set structure_has_changed=false to prevent another full LU if same nz_structure
		structure_has_changed = false;

	  }else if(values_have_changed || do_klu_refactor){
		clock_t start,finish;
		start = clock();

		klu_z_refactor ( Ap, Ai, reinterpret_cast <double*>(Ax), Symbolic, Numeric, &Common ) ;

	      	finish = clock();
		time_to_do_klu_refactor+= (double(finish)-double(start))/CLOCKS_PER_SEC;
		number_of_klu_refactor++;

		//we have done LU so set values_has_changed=false to prevent another LU if same values
		values_have_changed = false;
	  }
	  
	  	  
	  //Now do the F/B substitution
	  int result = klu_z_solve(Symbolic, Numeric, this->rows, Nrhs, (double*)(*RHS), &Common);	  
	  if(!result){
		std::cerr<<"Cannot do F/B substitution";
	  }
	  
      }
      
      
      
      SBase<std::complex<double> >* solve(const SBase<std::complex<double> > &B) const{
		throw std::runtime_error("Please, code me solve(const SBase<double> &B)");
      }
      
      Sparse<std::complex<double> > add(const Sparse<std::complex<double> > &B) const{
		if(this->rows!=B.rows || this->cols!=B.cols){
			throw std::runtime_error("Can not add two sparse matrices with different sizes");
    		}
	
    		const_cast<Sparse<std::complex<double> >&>(B).create_ccs();
		const_cast<Sparse<std::complex<double> >*>(this)->create_ccs();

    		Sparse<std::complex<double> > result(this->rows, this->cols);

    		for (int i = 0; i < this->cols; i++) {
      			int an = this->Ap[i];
      			int bn = B.Ap[i];

     			while (an < this->Ap[i+1] && bn < B.Ap[i+1]) {
        			if (this->Ai[an] == B.Ai[bn]) {
          				result.put(this->Ai[an],i, ((this->Ax[an]) + (B.Ax[bn])) );
					bn++;
				}else {
				        result.put(this->Ai[an],i, this->Ax[an] );
				}
				an++;
				
      			}
      			while (an < this->Ap[i+1]) {
        			result.put(this->Ai[an],i, (this->Ax[an]) );
        			an++;
      			}
      			while (bn < B.Ap[i+1]) {
        			result.put(B.Ai[bn],i, (B.Ax[bn]) );
        			bn++;
      			}
    		}
  
    		return result;
      }

      Sparse<std::complex<double> > subtract(const Sparse<std::complex<double> > &B) const{
		if(this->cols!=B.cols){
			throw std::runtime_error("Can not add two sparse matrices with different sizes");
    		}
	
    		const_cast<Sparse<std::complex<double> >&>(B).create_ccs();
		const_cast<Sparse<std::complex<double> >*>(this)->create_ccs();
		
    		Sparse<std::complex<double> > result(this->rows, this->cols);

    		for (int i = 0; i < this->rows; i++) {
      			int an = this->Ap[i];
      			int bn = B.Ap[i];

      			while (an < this->Ap[i+1] && bn < B.Ap[i+1]) {
        			if (this->Ai[an] == B.Ai[bn]) {
          				result.put(this->Ai[an],i, ((this->Ax[an]) - (B.Ax[bn])) );
					bn++;
				}else {
				        result.put(this->Ai[an],i, this->Ax[an] );
				}
				an++;
				
      			}
      			while (an < this->Ap[i+1]) {
        			result.put(this->Ai[an],i, (this->Ax[an]) );
        			an++;
      			}
      			while (bn < B.Ap[i+1]) {
        			result.put(B.Ai[bn],i, (B.Ax[bn]) );
        			bn++;
      			}
    		}
  
    		return result;
      }
      

      Sparse<std::complex<double> >& operator-= (const Sparse<std::complex<double> > &A){
	    if(this->rows!=A.rows || this->cols!=A.cols){
		throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    }
	
	    create_ccs();
            const_cast<Sparse<std::complex<double> >&>(A).create_ccs();


	    for (int i = 0; i < A.rows; i++) {
		  int an = Ap[i];
		  int bn = A.Ap[i];

		  while (an < Ap[i+1] && bn < A.Ap[i+1]) {
		      if (Ai[an] == A.Ai[bn]) {
			  put(i, Ai[an], (Ax[an]) - (A.Ax[bn]));
			  an++;
			  bn++;
		      }else {
			  put(i, A.Ai[bn], (A.Ax[bn]));
			  bn++;
		      }
		  }
	   
		  while (bn < A.Ap[i+1]) {
			put(i, A.Ai[bn], (A.Ax[bn]));
			bn++;
		  }
	    }
  
  	    //we have to use klu
	    return *this;
      }

       SBase<std::complex<double> >& operator -=(const SBase<std::complex<double> > &A){
  
		const Sparse<std::complex<double> >* BB = dynamic_cast<const Sparse<std::complex<double> >*>(&A);

  		if (BB) {
			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
    			return *this-=*BB;
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
	}
      
      Sparse<std::complex<double> >& operator+=(const Sparse<std::complex<double> > &A){
	    if(this->rows!=A.rows || this->cols!=A.cols){
		throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    }
	
	    create_ccs();
	    const_cast<Sparse<std::complex<double> >&>(A).create_ccs();


	    for (int i = 0; i < A.rows; i++) {
		  int an = Ap[i];
		  int bn = A.Ap[i];

		  while (an < Ap[i+1] && bn < A.Ap[i+1]) {
		      if (Ai[an] == A.Ai[bn]) {
			  put(i, Ai[an], ((Ax[an]) + A.Ax[bn]) );
			  an++;
			  bn++;
		      }else {
			  put(i, A.Ai[bn], A.Ax[bn] );
			  bn++;
		      }
		  }
	   
		  while (bn < A.Ap[i+1]) {
			put(i, A.Ai[bn], A.Ax[bn] );
			bn++;
		  }
	    }
  
	    return *this;
      }
     
      SBase<std::complex<double> >& operator +=(const SBase<std::complex<double> > &A){
  
		const Sparse<std::complex<double> >* BB = dynamic_cast<const Sparse<std::complex<double> >*>(&A);

  		if (BB) {
			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
    			return *this+=*BB;
  		}else{
    			throw std::runtime_error("I can only handle addition of Sparse+Sparse");
  		}
	}
      
      Sparse<std::complex<double> > operator+ (const Sparse<std::complex<double> > &B) const{
		return add(B);
      }

      SBase<std::complex<double> >* operator+(SBase<std::complex<double> > const &B) const{
  		const Sparse<std::complex<double> >* BB = dynamic_cast<const Sparse<std::complex<double> >*>(&B);

  		if (BB) {
    			return new Sparse<std::complex<double> >( *this + *BB );
  		}else{
    			throw std::runtime_error("I can only handle addition of Sparse+Sparse");
  		}
	}

      Sparse<std::complex<double> > operator- (const Sparse<std::complex<double> > &B) const{
		return subtract(B);
      }

      SBase<std::complex<double> >* operator-(SBase<std::complex<double> > const &B) const{
  		const Sparse<std::complex<double> >* BB = dynamic_cast<const Sparse<std::complex<double> >*>(&B);

  		if (BB) {
    			return new Sparse<std::complex<double> >( *this + *BB  );
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
      }
	  
      Sparse<std::complex<double> > operator* (const Sparse<std::complex<double> > &B) const{
		if (this->cols != B.rows) {
      			throw std::runtime_error("Can not multibly two sparse matrices with different rows and columns");
    		}

    		const_cast<Sparse<std::complex<double> >&>(*this).create_ccs();
    		const_cast<Sparse<std::complex<double> >&>(B).create_ccs();
    

    		int m = this->rows;
    		int n = B.cols;
 
    		Sparse<std::complex<double> > result(m, n);

    		for (int i = 0; i < B.cols; i++) {
      			for (int p = B.Ap[i]; p < B.Ap[i+1]; p++) {
        			int j = B.Ai[p];
        			for (int q = this->Ap[j]; q < this->Ap[j+1]; q++) {
          				result.add_to_entry(this->Ai[q], i, (this->Ax[q] * B.Ax[p]));
        			}
      			}
    		}

    		return result; 
      }

      SBase<std::complex<double> >* operator*(SBase<std::complex<double> > const &B) const{
  		const Sparse<std::complex<double> >* BB = dynamic_cast<const Sparse<std::complex<double> >*>(&B);

  		if (BB) {
    			return new Sparse<std::complex<double> >( *this * *BB  );
  		}else{
    			throw std::runtime_error("I can only handle multiblication of Sparse+Sparse");
  		}
	}

	Dense<std::complex<double> > operator * (Dense<std::complex<double> >& B){
	        if (this->cols != B.get_number_of_rows()) {
			throw std::runtime_error("Can not multible Sparse*Dense with inconsistence sizes");
		}

		create_ccs();
		
		Dense<std::complex<double> > result (this->rows, B.get_number_of_cols());

		
		for (int j = 0; j < B.get_number_of_cols(); j++) {	    
			    for (int i = 0; i < this->cols; i++) {
				  if (Ap[i] < Ap[i+1]) {
				  int an = this->Ap[i];
                    
				  while (an < this->Ap[i+1]) {
					result.add_to_entry (this->Ai[an] ,j, this->Ax[an] * B.get(i,j) );
					an++;
				  }
                    
				  
			    }
			}
		}

		return result;
	}
	
	
	Sparse<std::complex<double> >& operator *=(const Sparse<std::complex<double> > &A){
		
	
		return *this;
	}

       SBase<std::complex<double> >& operator *=(const SBase<std::complex<double> > &A){
		
  
		const Sparse<std::complex<double> >* BB = dynamic_cast<const Sparse<std::complex<double> >*>(&A);

  		if (BB) {
    			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
			return *this *= *BB;
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
	}
	
	Sparse<std::complex<double> > operator/(std::complex<double> val)const{
	    
	    Sparse<std::complex<double> > result = *this;
	  
	    result.create_ccs();
	    
	    for(int i=0; i<nnz; i++){
		result.Ax[i]/= val;
	    }
	    
	    return result;
	}

        Sparse<std::complex<double> > operator*(std::complex<double> val)const{
	    
	    Sparse<std::complex<double> > result = *this;
	  
	    result.create_ccs();
	    
	    for(int i=0; i<nnz; i++){
		result.Ax[i]*= val;
	    }
	    
	    return result;
	}
     
	  friend std::ostream& operator << (std::ostream &out , Sparse<std::complex<double> >& B);

	  Dense<std::complex<double> > operator + (BMatrix::Dense<std::complex<double> >& A){
		  if (A.get_number_of_rows() != this->rows
		      || A.get_number_of_cols() != this->cols) {
			    throw std::runtime_error("Can not Add Dense+Sparse with different sizes");
		  }

		  create_ccs();
		  
		  BMatrix::Dense<std::complex<double> > result = A;
    
    
		  for (int i = 0; i < this->cols; i++) {
		      int bn = Ap[i];

		      while (bn < Ap[i+1]) {
			    result.add_to_entry(Ai[bn] , i, Ax[bn] );
			    bn++;
		      }
		  }
		  
		  return result;
	  }

	  void report_timing(){
		std::cout<<"Number of klu_analyze done on this matrix = "<<number_of_klu_analyze <<std::endl;
		std::cout<<"Time to do all klu_analyze= "<<time_to_do_klu_analyze << std::endl;
      	        
		std::cout<<"Number of klu_fatcor done on this matrix = "<<number_of_klu_factor <<std::endl;
		std::cout<<"Time to do all klu_factor= "<<time_to_do_klu_factor <<std::endl;

		std::cout<<"Number of klu_refatcor done on this matrix = "<<number_of_klu_refactor <<std::endl;
		std::cout<<"Time to do all klu_refactor= "<<time_to_do_klu_refactor <<std::endl;
	  }

      protected:
      void runtime_error(const char* arg1);
};
  
  
inline std::ostream& operator << (std::ostream &out , Sparse<std::complex<double> >& B){

      const_cast< Sparse<std::complex<double> >&>(B).create_ccs();
      
      out<<"Ap=[";
      for(int i=0; i<B.cols+1; i++){
	  out<<B.Ap[i]<<" ";
      }
      out<<"]"<<std::endl;
      
      out<<"Ai=[";
      for(int i=0; i<B.nnz; i++){
	  out<<B.Ai[i]<<" ";
      }
      out<<"]"<<std::endl;
      
      out<<"Ax=[";
      for(int i=0; i<B.nnz; i++){
	  out<<(B.Ax[i])<<" ";
      }
      out<<"]"<<std::endl;
      
      return out;
} 

  

  
};

#endif