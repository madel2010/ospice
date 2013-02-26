
#ifndef BSPARSE_H
#define BSPARSE_H

#include <iostream>
#include <list>
#include <functional>
#include <algorithm>
#include <vector>

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "MWrap.h"
#include "klu.h"

#include <time.h>
#include <sys/time.h>
#include <sys/times.h>

namespace BMatrix{
  
template <class T>
class SBase{

protected:
	int rows; //number of rows
	int cols; //number of columns			
public:
	int get_number_of_rows(){ return rows;}
	int get_number_of_cols(){return cols;}

	virtual void put(int m, int n, T value)=0; 
	virtual T get(int m, int n) const=0;
};

template<class T>
class Sparse: public SBase<T>{

private:
      int nnz; //this will save the number of nonzeros
      bool ccs_created;
      
      //For the CCS 
      int* Ap;  //the pointer to columns position
      int* Ai; //the rows data
      T** Ax; 		//The values
      
      //For KLU
      BMatrix::MWrap<double>* MWrap_Ax; //The Mwrap class that holdes tha Ax values. It is used in the KLU routines
      klu_symbolic* Symbolic;
      klu_common Common;
      klu_B_numeric* Numeric;
      bool calculated_LU;
      bool calculated_Sparse_ordering;
      
      struct SparseElement{
	  int row;
	  T value;
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
	  Ax = new  T* [nnz]; //here we are using a pointer to T, because we want Ax to just point to the value in SparseElement.value instead of doing Ax[k] = SparseElement[k].value which will invoke operator = 
	  
	  
	  typename std::list<SparseElement>::iterator row_iter;
	  int k = 0;
	  

	  for(int i=0; i< this->cols; i++){
	      Ap[i+1] = Ap[i] + this->cols_lists[i].size();
	      
	      
	      for(row_iter = this->cols_lists[i].begin(); row_iter!=this->cols_lists[i].end(); row_iter++){
		  Ai[k] = row_iter->row;
		  Ax[k] = &row_iter->value; //note: Ax is a double pointer. that way we can only point to the value in row_iter instead of doing a copy constructor here
		  
		  k++;
	      }
	      
	      
	  }  

	  ccs_created = true;
      }
      
      void create_MWrap(){
	  
	  if(nnz==0) return;
	  
	  //First check that T is a block not scalar and it is Dense
	  if (!dynamic_cast<Dense<double>*>(Ax[0])) {
	      throw std::runtime_error("Solving a Block matrix with non-Dense blocks");
	  } 
	  
	  //This function should put the Ax vector in the Mwrap format to be ready for KLU routines
	  MWrap_Ax = new  BMatrix::MWrap<double>[nnz];
	  for(int i=0; i<nnz; i++){
	      MWrap_Ax[i] = dynamic_cast<DBase<double>*>(Ax[i]);
	  }
	  
      }
      
public:
      Sparse(){
	  this->rows=0; //Number of rows
	  this->cols=0; //Number of cols
	  
	  Ap = NULL;
	  Ai = NULL;
	  Ax = NULL;
	  MWrap_Ax = NULL;
	  
	  klu_defaults(&Common);
	  Common.scale=0;

	  ccs_created = false; 
	  calculated_LU = false;
	  calculated_Sparse_ordering = false;
	  
	  nnz = 0;
      }
      
      Sparse(int m, int n){
	  this->rows=m; //Number of rows
	  this->cols=n; //Number of cols
	  
	  cols_lists = new std::list<SparseElement>[n];
	  
	  Ap = new int[n+1]; //We know that Ap is always (columns_no+1)
	  Ai = NULL; //we still do not know the rows
	  Ax = NULL; //we still do not know the values
	  MWrap_Ax = NULL;
	  klu_defaults(&Common);
	  Common.scale=0;

	  ccs_created = false;
	  calculated_LU = false;
	  calculated_Sparse_ordering = false;
	  
	  nnz = 0;
      }
      
     
      
      //copy constructor
      Sparse(const Sparse<T>&A){
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
	  MWrap_Ax = NULL;
	  calculated_LU = false;
	  calculated_Sparse_ordering = false;
	  
	  klu_defaults(&Common);
	  Common.scale=0;

      }
      
	Sparse<T>* clone(){
		return new Sparse<T>(*this);
	}
	
	void create(int m, int n){ 
		this->rows=m; //Number of rows
		this->cols=n; //Number of cols
	  
		cols_lists = new std::list<SparseElement>[n];
	  
		Ap = new int[n+1]; //We know that Ap is always (columns_no+1)
		Ai = NULL; //we still do not know the rows
		Ax = NULL; //we still do not know the values
		MWrap_Ax = NULL;
		klu_defaults(&Common);
	  	Common.scale=0;

		ccs_created = false;
		calculated_LU = false;
		calculated_Sparse_ordering = false;
	  
		nnz = 0;
	}
	
      ~Sparse(){
	delete[] Ap;
	Ap=NULL;
	delete[] Ai;
	Ai = NULL;	
      }

      Sparse<T>& operator=(T val){}
      
      SBase<T>* scale(T* val)const {}
	
      Sparse<T>& operator/=(T val){throw std::runtime_error("Please code me Sparse::operator/= (T val)");}
      
      
      
      Sparse<T>& operator=(const Sparse<T>&A){
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
	  calculated_LU = false;
	  calculated_Sparse_ordering = false;
	  
	 return *this;
      }
      
      //Get a value at row m and column n
      T get(int m, int n) const{
	  T result;
	  
	  //search column for the row
	  typename std::list<SparseElement>::iterator row_iterator;
	  
	 row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	 
	 if(row_iterator!=cols_lists[n].end()){ //we have found the row
	    result = row_iterator->value;
	    
	 }else{
	    result = 0;
	 }
	 
	 return result;
	  
      }
      
      double det(){

      }

      bool is_zero(){

      }

      //add a  value in row m and column n to the existing value
      void add(int m, int n, T value){
	  
	  //search column for the row in case we have added the row already
	  typename std::list<SparseElement>::iterator row_iterator;
	  
	 row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	 
	 
	 if(row_iterator!=cols_lists[n].end()){ //we have found a value already
	    row_iterator->value += value;
	    
	 }else{  //we have no value at this row. We have to add the row first
	      row_iterator = cols_lists[n].begin();
	      while(true){
		  if(row_iterator->row > m){  //we have found the first row greater than the required row
		    
			 SparseElement new_element; //create a new element structure
			 new_element.row = m;  //add the row to the new element structure
			 new_element.value = value; //add the value to the new element structure
			 
			 cols_lists[n].insert(row_iterator , new_element); 
			 
			 nnz++; //increase the number of non zeros in the matrix
			 
			 break;
		  }
		  row_iterator++;
	      }
	      ccs_created = false;
	}
	 
	 
	 calculated_LU = false;
	 calculated_Sparse_ordering = false;
	 
      }
      
      //put value in row m and column n
      void put(int m, int n, T value){
	  
	  //search column for the row in case we have added the row already
	  typename std::list<SparseElement>::iterator row_iterator;
	  
	 row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	 
	 
	 if(row_iterator!=cols_lists[n].end()){ //we have found a value already
	    row_iterator->value = value;
	    
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
	      ccs_created = false;
	 }
	 
	 
	 calculated_LU = false;
	 calculated_Sparse_ordering = false;
	 
      }

      
      //get the number of non zeros entries
      int get_nnz(){ return nnz; }
      
      void solve(BMatrix::Dense<T>& RHS, const int Nrhs=1){
	  
	  create_ccs();
	  create_MWrap();
	  if(!calculated_Sparse_ordering){
		Symbolic = klu_analyze(this->rows, Ap, Ai, &Common) ;
	  }
	  
	  //change RHS to MWrap
	  BMatrix::MWrap<double>* MWrap_RHS = new  BMatrix::MWrap<double>[RHS.get_number_of_rows()];
	  for(int i=0; i<RHS.get_number_of_rows(); i++){
	      MWrap_RHS[i] = RHS.get_pointer(i,0);
	  }


	  //Do LU factorization if not done before
	  if(!calculated_LU){
	      Numeric  = klu_B_factor ( Ap, Ai, MWrap_Ax, Symbolic, &Common ) ;
	  }
	  
	  
	  
	  //Now do the F/B substitution
	  int result = klu_B_solve(Symbolic, Numeric, this->rows, Nrhs, MWrap_RHS, &Common);	  
	  if(!result){
		std::cerr<<"Cannot do F/B substitution";
	  }

	  sparse_free_numeric();
      }
      
      void sparse_free_numeric(){
	  klu_B_free_numeric(&Numeric ,  &Common);
      }
      
      SBase<T>* solve(const SBase<T> &B) const{
		throw std::runtime_error("Please, code me solve(const SBase<T> &B)");
      }
      
      Sparse<T> add(const Sparse<T> &B) const{
		if(this->rows!=B.rows || this->cols!=B.cols){
			throw std::runtime_error("Can not add two sparse matrices with different sizes");
    		}
	
    		const_cast<Sparse<T>&>(B).create_ccs();
		const_cast<Sparse<T>*>(this)->create_ccs();
		
    		Sparse<T> result(this->rows, this->cols);

    		for (int i = 0; i < this->cols; i++) {
      			int an = this->Ap[i];
      			int bn = B.Ap[i];

     			while (an < this->Ap[i+1] && bn < B.Ap[i+1]) {
        			if (this->Ai[an] == B.Ai[bn]) {
          				result.put(this->Ai[an],i, (*(this->Ax[an]) + *(B.Ax[bn])) );
					bn++;
				}else  {
				        result.put(this->Ai[an],i, *(this->Ax[an]) );
				}
				an++;
				
      			}
      			while (an < this->Ap[i+1]) {
        			result.put(this->Ai[an],i, *(this->Ax[an]) );
        			an++;
      			}
      			while (bn < B.Ap[i+1]) {
        			result.put(B.Ai[bn],i, *(B.Ax[bn]) );
        			bn++;
      			}
    		}
  
    		return result;
      }

      Sparse<T> subtract(const Sparse<T> &B) const{
		if(this->cols!=B.cols){
			throw std::runtime_error("Can not add two sparse matrices with different sizes");
    		}
	
    		const_cast<Sparse<T>&>(B).create_ccs();
		const_cast<Sparse<T>*>(this)->create_ccs();
		
    		Sparse<T> result(this->rows, this->cols);

    		for (int i = 0; i < this->rows; i++) {
      			int an = this->Ap[i];
      			int bn = B.Ap[i];

      			while (an < this->Ap[i+1] && bn < B.Ap[i+1]) {
        			if (this->Ai[an] == B.Ai[bn]) {
          				result.put(this->Ai[an],i, (*(this->Ax[an]) + *(B.Ax[bn])) );
					bn++;
				}else  {
				        result.put(this->Ai[an],i, *(this->Ax[an]) );
				}
				an++;
				
      			}
      			while (an < this->Ap[i+1]) {
        			result.put(this->Ai[an],i, *(this->Ax[an]) );
        			an++;
      			}
      			while (bn < B.Ap[i+1]) {
        			result.put(B.Ai[bn],i, *(B.Ax[bn]) );
        			bn++;
      			}
    		}
  
    		return result;
      }
      

      Sparse<T>& operator-= (const Sparse<T> &A){
	    if(this->rows!=A.rows || this->cols!=A.cols){
		throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    }
	
	    create_ccs();
            const_cast<Sparse<T>&>(A).create_ccs();


	    for (int i = 0; i < A.rows; i++) {
		  int an = Ap[i];
		  int bn = A.Ap[i];

		  while (an < Ap[i+1] && bn < A.Ap[i+1]) {
		      if (Ai[an] == A.Ai[bn]) {
			  put(i, Ai[an], (*Ax[an]) - (*A.Ax[bn]));
			  an++;
			  bn++;
		      }else {
			  put(i, A.Ai[bn], (*A.Ax[bn]));
			  bn++;
		      }
		  }
	   
		  while (bn < A.Ap[i+1]) {
			put(i, A.Ai[bn], (*A.Ax[bn]));
			bn++;
		  }
	    }
  
	    return *this;
      }

       SBase<T>& operator -=(const SBase<T> &A){
  
		const Sparse<T>* BB = dynamic_cast<const Sparse<T>*>(&A);

  		if (BB) {
			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
    			return *this-=*BB;
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
	}
      
      Sparse<T>& operator+=(const Sparse<T> &A){
	    if(this->rows!=A.rows || this->cols!=A.cols){
		throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    }
	
	    create_ccs();
	    const_cast<Sparse<T>&>(A).create_ccs();


	    for (int i = 0; i < A.rows; i++) {
		  int an = Ap[i];
		  int bn = A.Ap[i];

		  while (an < Ap[i+1] && bn < A.Ap[i+1]) {
		      if (Ai[an] == A.Ai[bn]) {
			  put(i, Ai[an], ((*Ax[an]) + (*A.Ax[bn])) );
			  an++;
			  bn++;
		      }else {
			  put(i, A.Ai[bn], (*A.Ax[bn]) );
			  bn++;
		      }
		  }
	   
		  while (bn < A.Ap[i+1]) {
			put(i, A.Ai[bn], (*A.Ax[bn]) );
			bn++;
		  }
	    }
  
	    return *this;
      }
     
      SBase<T>& operator +=(const SBase<T> &A){
  
		const Sparse<T>* BB = dynamic_cast<const Sparse<T>*>(&A);

  		if (BB) {
			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
    			return *this+=*BB;
  		}else{
    			throw std::runtime_error("I can only handle addition of Sparse+Sparse");
  		}
	}
      
      Sparse<T> operator+ (const Sparse<T> &B) const{
		return add(B);
      }

      SBase<T>* operator+(SBase<T> const &B) const{
  		const Sparse<T>* BB = dynamic_cast<const Sparse<T>*>(&B);

  		if (BB) {
    			return new Sparse<T>( *this + *BB );
  		}else{
    			throw std::runtime_error("I can only handle addition of Sparse+Sparse");
  		}
	}

      Sparse<T> operator- (const Sparse<T> &B) const{
		return subtract(B);
      }

      SBase<T>* operator-(SBase<T> const &B) const{
  		const Sparse<T>* BB = dynamic_cast<const Sparse<T>*>(&B);

  		if (BB) {
    			return new Sparse<T>( *this + *BB  );
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
      }
	  
      Sparse<T> operator* (const Sparse<T> &B) const{
		if (this->cols != B.rows) {
      			throw std::runtime_error("Can not multibly two sparse matrices with different rows and columns");
    		}

    		const_cast<Sparse<T>&>(*this).create_ccs();
    		const_cast<Sparse<T>&>(B).create_ccs();
    

    		int m = this->rows;
    		int n = B.cols;
 
    		Sparse<T> result(m, n);

    		for (int i = 0; i < B.cols; i++) {
      			for (int p = B.Ap[i]; p < B.Ap[i+1]; p++) {
        			int j = B.Ai[p];
        			for (int q = this->Ap[j]; q < this->Ap[j+1]; q++) {
          				result.add(this->Ai[q], i, ((*this->Ax[q])*(*B.Ax[p])));
        			}
      			}
    		}

    		return result; 
      }

      SBase<T>* operator*(SBase<T> const &B) const{
  		const Sparse<T>* BB = dynamic_cast<const Sparse<T>*>(&B);

  		if (BB) {
    			return new Sparse<T>( *this * *BB  );
  		}else{
    			throw std::runtime_error("I can only handle multiblication of Sparse+Sparse");
  		}
	}

	Dense<T> operator * (Dense<T>& B){
	        if (this->cols != B.get_number_of_rows()) {
			throw std::runtime_error("Can not multible Sparse*Dense with inconsistence sizes");
		}

         
		Dense<T> result (this->rows, B.get_number_of_cols());

		
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
	
	
	Sparse<T>& operator *=(const Sparse<T> &A){
		
	
		return *this;
	}

       SBase<T>& operator *=(const SBase<T> &A){
		
  
		const Sparse<T>* BB = dynamic_cast<const Sparse<T>*>(&A);

  		if (BB) {
    			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
			return *this *= *BB;
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
	}
	     
	Sparse<T> operator/(T val){
		create_ccs();
		
		for(int i=0; i<nnz; i++){
		    (*Ax[i]) /= val;
		}
	}
	     
	template<class U>
	      friend std::ostream& operator << (std::ostream &out , Sparse<U>& B);
	
	Dense<T> operator + (BMatrix::Dense<T>& A ){
	      if (A.get_number_of_rows() != this->rows
		      || A.get_number_of_cols() != this->cols) {
		  throw std::runtime_error("Can not Add Dense+Sparse with different sizes");
	      }

	      BMatrix::Dense<T> result = A;
    
    
	      for (int i = 0; i < this->cols; i++) {
		    int bn = Ap[i];

	      while (bn < Ap[i+1]) {
		    result.add_to_entry(Ai[bn] , i, (*Ax[bn]) );
		    bn++;
	      }
	    }
	    
	    return result;
	}
	      
	      
      protected:
      void runtime_error(const char* arg1);
};
  
  
template<class U> std::ostream& operator << (std::ostream &out , Sparse<U>& B){

      const_cast<Sparse<U>&>(B).create_ccs();
      
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
	  out<<(*B.Ax[i])<<" ";
      }
      out<<"]"<<std::endl;
      
      return out;
} 

  

  
};

#endif