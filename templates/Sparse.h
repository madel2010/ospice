
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
#include "DoubleDense.h"

#include <time.h>
#include <sys/time.h>
#include <sys/times.h>

namespace BMatrix{
  
template <class T>
class SBase{

      
protected:
	int rows; //number of rows
	int cols; //number of columns	
	
	T zero_element; //we want to know what a zero element is, this is usefull for the case of matrices of matrices when we use put anf get
	
public:
     
	int get_number_of_rows()const{ return rows;}
	int get_number_of_cols()const{return cols;}

	virtual void put(int m, int n, T value)=0; 
	virtual T get(int m, int n) const=0;

	SBase(T _zero_element){zero_element=_zero_element;}
	SBase(){zero_element=0.0;};
	
	virtual ~SBase(){}

};

template<class T>
class Sparse: public SBase<T>{

private:
          
      int nnz; //this will save the number of nonzeros
      bool ccs_created;
      bool matrix_created;

      //For the CCS 
      int* Ap;  //the pointer to columns position
      int* Ai; //the rows data
      T* Ax; 		//The values
      
      
      //For KLU
      klu_symbolic* Symbolic;
      klu_common Common;
      void* Numeric;
      
      bool structure_has_changed;
      bool values_have_changed;
	
      BMatrix::MWrap<double>* MWrap_Ax; //The Mwrap class that holdes tha Ax values. It is used in the KLU routines
       
      //the following numbers are for storing the time for doing sparse order and LU factorization.
      double time_to_do_klu_analyze;
      double time_to_do_klu_factor;
      double time_to_do_klu_refactor;
      int number_of_klu_analyze;
      int number_of_klu_factor;
      int number_of_klu_refactor;

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
      int* first_row; //save the first row we already added for each column. Very useful to speedup put()
      int* last_row; //save the last row we already added for each column. Very useful to speedup put()
      typename std::list<SparseElement>::iterator *Last_accessed_ele_in_col; //This keeps track of the last element accessed in each col. This is useful becuase for example in the multiplicatio we are moving element by element, so we know that the element to be currently accessed is very close to the last element.
            
      //Create the CCS formal from the linked list
      void create_ccs(){
	  if(this->ccs_created==true){
		return;  //we do not need to redo it if we already have found the css
	  }  
	  
	  this->Ap[0] = 0; //this->Ap[0] is always zero
	  
	  this->Ai = new  int[this->nnz];
	  this->Ax = new  T[this->nnz]; 
	  
	  
	  //typedef typename std::list<SparseElement>::iterator row_iter;
	  typename std::list<SparseElement>::iterator row_iter;
	  int k = 0;
	  

	  for(int i=0; i< this->cols; i++){
	      this->Ap[i+1] = this->Ap[i] + this->cols_lists[i].size();
	      
	      
	      for(row_iter = this->cols_lists[i].begin(); row_iter!=this->cols_lists[i].end(); row_iter++){
		  this->Ai[k] = row_iter->row;
		  this->Ax[k] = row_iter->value; 
		  
		  k++;
	      }
	      
	      
	  }  
	  
	  this->ccs_created = true;
	  
      }
      
  
      //search col if a row has been added already
      typename std::list<SparseElement>::iterator search_column(int col, int row){
	  typename std::list<SparseElement>::iterator result;
	
	  //we have to get the last accessed element in the col
	  if(row >= (*Last_accessed_ele_in_col[col]).row){ //this means that we have to search forward
		result = std::find_if( Last_accessed_ele_in_col[col], cols_lists[col].end(), std::bind2nd( FindRow(), row ) );
	  }else if(row < (*Last_accessed_ele_in_col[col]).row){ //we have to search backward

		typename std::list<SparseElement>::reverse_iterator search_results;

		typename std::list<SparseElement>::reverse_iterator start(Last_accessed_ele_in_col[col]);
		search_results = std::find_if( start , cols_lists[col].rend(), std::bind2nd( FindRow(), row ) );

		result = --search_results.base();
	  }	

	  return result;


      }

      void create_MWrap(){
	  
	  if(this->nnz==0) return;
	  
	  //First check that T is a block not scalar and it is Dense
	  if (!dynamic_cast<Dense<double>*>(&(this->Ax[0]))) {
	      throw std::runtime_error("Solving a Block matrix with non-Dense blocks");
	  } 
	  
	  //This function should put the Ax vector in the Mwrap format to be ready for KLU routines
	  MWrap_Ax = new  BMatrix::MWrap<double>[this->nnz];
	  for(int i=0; i<this->nnz; i++){
	      MWrap_Ax[i] = dynamic_cast<DBase<double>*>(&(this->Ax[i]));
	  }
      }
      
public:
      Sparse(){
	  this->rows=0; //Number of rows
	  this->cols=0; //Number of cols
	  this->nnz = 0;

	  this->Ap = NULL;
	  this->Ai = NULL;
	  this->Ax = NULL;
	  this->Symbolic = NULL;
	  this->Numeric = NULL;
	  

	  cols_lists = NULL;
	  first_row = NULL;
	  last_row = NULL;
	  Last_accessed_ele_in_col = NULL;
 
	  klu_defaults(&this->Common);
	  this->Common.scale=0;

	  this->matrix_created = false;
	  this->ccs_created = false; 
	  this->structure_has_changed = true;
	  this->values_have_changed = true;

	  

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
      }
      
      Sparse(int m, int n){
	  /*this->rows=m; //Number of rows
	  this->cols=n; //Number of cols
	  this->nnz = 0;

	  cols_lists = new std::list<SparseElement>[n];
	  first_row = new int[n];
	  last_row = new int[n];
	  memset(first_row, -1 , n*sizeof(int));
	  memset(last_row, -1 , n*sizeof(int));
	  Last_accessed_ele_in_col = new typename std::list<SparseElement>::iterator[n];
	  

	  this->Ap = new int[n+1]; //We know that this->Ap is always (columns_no+1)
	  this->Ai = NULL; //we still do not know the rows
	  this->Ax = NULL; //we still do not know the values
	  this->Symbolic = NULL;
	  this->Numeric = NULL;
	  
	  klu_defaults(&this->Common);
	  this->Common.scale=0;

	  this->matrix_created = true;
	  this->ccs_created = false;
	  this->structure_has_changed = true;
	  this->values_have_changed = true;

	  

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
	  */
	  this->matrix_created = false;
	  create(m, n);
      }
      
      //copy constructor
      Sparse(const Sparse<T >&A){
	  
	  int n = A.cols;	

	  cols_lists = new std::list<SparseElement>[n];
	  first_row = new int[n];
	  last_row = new int[n];
	  
	  memcpy(this->first_row , A.first_row , n*sizeof(int));
	  memcpy(this->last_row , A.last_row , n*sizeof(int));
	  
	  Last_accessed_ele_in_col = new typename std::list<SparseElement>::iterator[n];

	  this->Ap = new int[n+1]; //We know that this->Ap is always (columns_no+1)
	  this->Ai = NULL; //we still do not know the rows
	  this->Ax = NULL; //we still do not know the values
	  this->Symbolic = NULL;
          this->Numeric = NULL;
	  

	  this->nnz = A.nnz;
	  this->rows=A.rows; //Number of rows
	  this->cols=A.cols; //Number of cols
	  
	  
	  

	  for(int i=0 ; i< A.cols; i++){
	      cols_lists[i] = A.cols_lists[i];
	      Last_accessed_ele_in_col[i] =  cols_lists[i].begin();
	  }
	  
	  
	  
	  //TODO copy the CCS structure as well
	  this->ccs_created = false;
	
	  this->matrix_created = true;

	  this->structure_has_changed = true;
	  this->values_have_changed = true;

	  this->Symbolic = NULL;
	  this->Numeric = NULL;

	  klu_defaults(&this->Common);
	  this->Common.scale=0;

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
      }

      Sparse(int m, int n,T zero_element):SBase<T>(zero_element){
	  this->matrix_created = false;
	  create(m, n);
      }

       //convert dense to sparse
       Sparse(const Dense<T >&A){
	  int n = A.get_number_of_cols();	
	  int m = A.get_number_of_rows();

	  cols_lists = new std::list<SparseElement>[n];
	  first_row = new int[n];
	  last_row = new int[n];
	  
	  memset(first_row, -1 , n*sizeof(int));
	  memset(last_row, -1 , n*sizeof(int));
	  
	  Last_accessed_ele_in_col = new typename std::list<SparseElement>::iterator[n];
	  	  

	  this->Ap = new int[n+1]; //We know that this->Ap is always (columns_no+1)
	  this->Ai = NULL; //we still do not know the rows
	  this->Ax = NULL; //we still do not know the values
	  this->Symbolic = NULL;
          this->Numeric = NULL;
	  this->matrix_created = true;
	  
	  this->nnz = 0;

          for(int i=0 ; i< A.get_number_of_cols(); i++){
	      for (int j=0;j<A.get_number_of_rows(); j++){
			T value = A.get(j,i);
			if(value!=0.0){
	        		SparseElement new_element; //create a new element structure
	  			new_element.row = j;  //add the row to the new element structure
	  			new_element.value = value;  //add the value to the new element structure
		 		cols_lists[i].push_back(new_element); 	

					
				if(j> last_row[i]){
					last_row[i] = j;
				}else if(j < first_row[i]){
					first_row[i] = j;
				}
	
				this->nnz++;
			}
	      }
	      Last_accessed_ele_in_col[i] = cols_lists[i].begin();
	  }

	    
	  
	  //TODO copy the CCS structure as well
	  this->ccs_created = false;
	
	  this->structure_has_changed = true;
	  this->values_have_changed = true;

	  klu_defaults(&this->Common);
	  this->Common.scale=0;

	  //the time report data
	  time_to_do_klu_analyze = 0;
      	  time_to_do_klu_factor = 0;
	  time_to_do_klu_refactor = 0;
          number_of_klu_analyze = 0;
          number_of_klu_factor = 0;
	  number_of_klu_refactor = 0;
        }

	Sparse<T >* clone(){
		return new Sparse<T >(*this);
	}
	
	void create(int m, int n){ 

		if(this->matrix_created){
			return;
		}

		this->rows=m; //Number of rows
		this->cols=n; //Number of cols
	  
		cols_lists = new std::list<SparseElement>[n];
	  	first_row = new int[n];
	  	last_row = new int[n];
	  	memset(first_row, -1 , n*sizeof(int));
	  	memset(last_row, -1 , n*sizeof(int));
	        Last_accessed_ele_in_col = new typename std::list<SparseElement>::iterator[n];
	  

		this->Ap = new int[n+1]; //We know that this->Ap is always (columns_no+1)
		this->Ai = NULL; //we still do not know the rows
		this->Ax = NULL; //we still do not know the values
		this->Symbolic = NULL;
        	this->Numeric = NULL;
		this->matrix_created = true;

		klu_defaults(&this->Common);
	  	this->Common.scale=0;

		this->ccs_created = false;
		this->structure_has_changed = true;
	  	this->values_have_changed = true;
	  	
		this->nnz = 0;

		//the time report data
	  	time_to_do_klu_analyze = 0;
      	  	time_to_do_klu_factor = 0;
	  	time_to_do_klu_refactor = 0;
          	number_of_klu_analyze = 0;
          	number_of_klu_factor = 0;
	  	number_of_klu_refactor = 0;
	}
	
	void reset(){
	      typename std::list<SparseElement>::iterator row_iterator;

	      for(int col=0; col< this->cols; col++){
	    	  for(row_iterator= this->cols_lists[col].begin() ; row_iterator!=this->cols_lists[col].end(); row_iterator++){
			row_iterator->value = 0.0;
	    	  }
	      }
	      this->values_have_changed = true;
	}

      ~Sparse(){
	if(this->Ap) delete[] this->Ap;
	this->Ap=NULL;
	if(this->Ai) delete[] this->Ai;
	this->Ai = NULL;
	if(this->Ax) delete[] this->Ax;	

	if(cols_lists)	delete[] cols_lists;
	cols_lists = NULL;

	if(first_row) delete[] first_row;
	first_row = NULL;
        
	if(last_row) delete[] last_row;
	last_row = NULL;

        if(Last_accessed_ele_in_col) delete[] Last_accessed_ele_in_col;
	Last_accessed_ele_in_col = NULL;

	klu_free_symbolic (&this->Symbolic, &this->Common);
	this->Symbolic=NULL;
	
	klu_B_numeric* my_numeric = (klu_B_numeric*)this->Numeric;
        klu_B_free_numeric (&my_numeric, &this->Common); 
	
	this->Numeric=NULL;

      }

      Sparse<T >& operator=(T val){throw std::runtime_error("operator=(T val) not codded yet");}
      
      SBase<T >* scale(T* val)const {throw std::runtime_error("scale(T* val) not coded yet");}
	
      
      //operator =
      Sparse<T >& operator=(const Sparse<T >&A){
	  
	  if(this->rows!= A.rows || this->cols!= A.cols){
		throw std::runtime_error("Can not equate two sparse matrices with different size");
 	  }

	  int n = A.cols;
	  int m = A.rows;
	  if(!this->matrix_created){
	  	

	  	cols_lists = new std::list<SparseElement>[n];
	  	first_row = new int[n];
       	  	last_row = new int[n];
	  	
		Last_accessed_ele_in_col = new typename std::list<SparseElement>::iterator[n];

	  	this->Ap = new int[n+1]; //We know that this->Ap is always (columns_no+1)
	  	this->Ai = NULL; //we still do not know the rows
	  	this->Ax = NULL; //we still do not know the values
	  	this->Symbolic = NULL;
          	this->Numeric = NULL;
	  	this->matrix_created = true;
	  }

	  memcpy(this->first_row , A.first_row , n*sizeof(int));
	  memcpy(this->last_row , A.last_row , n*sizeof(int));
	  
	  this->nnz = A.nnz;
	  this->rows=A.rows; //Number of rows
	  this->cols=A.cols; //Number of cols

	  for(int i=0 ; i< A.cols; i++){
	      cols_lists[i] = A.cols_lists[i];
	      Last_accessed_ele_in_col[i] = cols_lists[i].begin();
	  }
	  
	  
	  //TODO copy the CCS structure as well
	  this->ccs_created = false; //for now we just make the new object create the CCS
	  this->structure_has_changed = true;
	  this->values_have_changed = true;

	 return *this;
      }
      
      //Get a value at row m and column n
      T get(int m, int n) const{
	  
	  
	  
	  //check if the row of the new value is greater/less than the last/first row we have already added
	  if(m > last_row[n] || m < first_row[n] || this->nnz==0){
	     return this->zero_element;
	  }
	  
	  T result;
	  
	  //search column for the row
	  //typename std::list<SparseElement>::iterator row_iterator;
	  typename std::list<SparseElement>::iterator row_iterator;

	 row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	 
	 if(row_iterator!=cols_lists[n].end()){ //we have found the row
	    result = row_iterator->value;
	    
	 }else{
	    result = this->zero_element;
	 }
	 
	 return result;
	  
      }
      
      double det(){
	      throw std::runtime_error("det() not coded yet");
      }

      bool is_zero(){
		if (this->nnz==0){
			return true;
		}else{
			return false;
		}	
      }

      //add a  value in row m and column n to the existing value
      void add_to_entry(int m, int n, T value){
	  
	   if(value==0.0) return;
	  
	  //check if this is the first element
	  if(this->nnz==0){
	        SparseElement new_element; //create a new element structure
		new_element.row = m;  //add the row to the new element structure
		new_element.value = value; //add the value to the new element structure
			 
		cols_lists[n].push_back(new_element); 	 
		this->nnz++; //increase the number of non zeros in the matrix	

		last_row[n] = m;
		first_row[n] = m;
		
		Last_accessed_ele_in_col[n] = --cols_lists[n].end() ;

		this->structure_has_changed = true;
		this->ccs_created = false;
	
		return ;
	  //check if the row of the new value is greater than the last row we have already added	
	  }else if(m > last_row[n]){
		SparseElement new_element; //create a new element structure
		new_element.row = m;  //add the row to the new element structure
		new_element.value = value; //add the value to the new element structure
			 
		cols_lists[n].push_back(new_element); 	 
		this->nnz++; //increase the number of non zeros in the matrix	

		last_row[n] = m;

		Last_accessed_ele_in_col[n] = --cols_lists[n].end() ;

		this->structure_has_changed = true;
		this->ccs_created = false;
	
		return ;
	  //check if the row of the new value is less than the last row we have already added
	  }else if(m < first_row[n]){ 
		SparseElement new_element; //create a new element structure
		new_element.row = m;  //add the row to the new element structure
		new_element.value = value; //add the value to the new element structure
			 
		cols_lists[n].push_front(new_element); 	 
		this->nnz++; //increase the number of non zeros in the matrix	

		first_row[n] = m;

		Last_accessed_ele_in_col[n] = cols_lists[n].begin();

		this->structure_has_changed = true;
		this->ccs_created = false;

		return ;
	  }

	  //search column for the row in case we have added the row already
	  //typename std::list<SparseElement>::iterator row_iterator;
	  typename std::list<SparseElement>::iterator row_iterator;

	  //row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	  row_iterator = search_column(n,m);
	 
	  if(row_iterator!=cols_lists[n].end()){ //we have found a value already
	     row_iterator->value += value;
	     this->values_have_changed = true;

	     Last_accessed_ele_in_col[n] = row_iterator;

	  }else{  //we have no value at this row. We have to add the row first
	       row_iterator = cols_lists[n].begin();
	      
	       bool found_larger_row = false;
	       while(row_iterator!=cols_lists[n].end()){
		       if(row_iterator->row > m){  //we have found the first row greater than the required row
		    
			     SparseElement new_element; //create a new element structure
			     new_element.row = m;  //add the row to the new element structure
			     new_element.value = value; //add the value to the new element structure
			 
			     cols_lists[n].insert(row_iterator , new_element); 
			 
			     this->nnz++; //increase the number of non zeros in the matrix
			 
			     Last_accessed_ele_in_col[n] = row_iterator;

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
			 
			  Last_accessed_ele_in_col[n] = --cols_lists[n].end();

			  this->nnz++; //increase the number of non zeros in the matrix
		}
	  
	
	       this->structure_has_changed = true;
	}
	  
	  if(this->ccs_created){
		if(this->Ai) delete[] this->Ai;
		this->Ai = NULL;
		if(this->Ax) delete[] this->Ax;
		this->Ax = NULL;
	  }
	  this->ccs_created = false; 
	 
      }
      
      //put value in row m and column n
      void put(int m, int n, T value){
	  
	  if(value==0.0) return;
	  
	  //check if this is the first element
	  if(this->nnz==0){
	        SparseElement new_element; //create a new element structure
		new_element.row = m;  //add the row to the new element structure
		new_element.value = value; //add the value to the new element structure
			 
		cols_lists[n].push_back(new_element); 	 
		this->nnz++; //increase the number of non zeros in the matrix	

		last_row[n] = m;
		first_row[n] = m;
		
		Last_accessed_ele_in_col[n] = --cols_lists[n].end() ;

		this->structure_has_changed = true;
		this->ccs_created = false;
	
		return ;
	  //check if the row of the new value is greater than the last row we have already added	
	  }else if(m > last_row[n]){
		SparseElement new_element; //create a new element structure
		new_element.row = m;  //add the row to the new element structure
		new_element.value = value; //add the value to the new element structure
			 
		cols_lists[n].push_back(new_element); 	 
		this->nnz++; //increase the number of non zeros in the matrix	

		last_row[n] = m;

		Last_accessed_ele_in_col[n] = --cols_lists[n].end() ;

		this->structure_has_changed = true;
		this->ccs_created = false;
	
		return ;
	  //check if the row of the new value is less than the last row we have already added
	  }else if(m < first_row[n]){ 
		SparseElement new_element; //create a new element structure
		new_element.row = m;  //add the row to the new element structure
		new_element.value = value; //add the value to the new element structure
			 
		cols_lists[n].push_front(new_element); 	 
		this->nnz++; //increase the number of non zeros in the matrix	

		first_row[n] = m;

		Last_accessed_ele_in_col[n] = cols_lists[n].begin();

		this->structure_has_changed = true;
		this->ccs_created = false;

		return ;
	  }

	  //search column for the row in case we have added the row already
	  //typename std::list<SparseElement>::iterator row_iterator;
	  typename std::list<SparseElement>::iterator row_iterator;

	  //row_iterator = find_if( cols_lists[n].begin(), cols_lists[n].end(), std::bind2nd( FindRow(), m ) );
	  row_iterator = search_column(n,m);
	 
	  if(row_iterator!=cols_lists[n].end()){ //we have found a value already
	     row_iterator->value = value;
	     this->values_have_changed = true;

	     Last_accessed_ele_in_col[n] = row_iterator;

	  }else{  //we have no value at this row. We have to add the row first
	       row_iterator = cols_lists[n].begin();
	      
	       bool found_larger_row = false;
	       while(row_iterator!=cols_lists[n].end()){
		       if(row_iterator->row > m){  //we have found the first row greater than the required row
		    
			     SparseElement new_element; //create a new element structure
			     new_element.row = m;  //add the row to the new element structure
			     new_element.value = value; //add the value to the new element structure
			 
			     cols_lists[n].insert(row_iterator , new_element); 
			 
			     this->nnz++; //increase the number of non zeros in the matrix
			 
			     Last_accessed_ele_in_col[n] = row_iterator;

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
			 
			  Last_accessed_ele_in_col[n] = --cols_lists[n].end();

			  this->nnz++; //increase the number of non zeros in the matrix
		}
	  
	
	       this->structure_has_changed = true;
	}
	  
	  if(this->ccs_created){
		if(this->Ai) delete[] this->Ai;
		this->Ai = NULL;
		if(this->Ax) delete[] this->Ax;
		this->Ax = NULL;
	  }
	  this->ccs_created = false; 
      }


 
      //get the number of non zeros entries
      int get_nnz(){ return this->nnz; }
      
      //This re writes the RHS
       void solve(BMatrix::Dense<T>& RHS, bool do_klu_refactor=false){
	  
	  create_ccs();
	  int Nrhs = RHS.get_number_of_cols();
	  
	  create_MWrap();
	  if(this->structure_has_changed){
		this->Symbolic = klu_analyze(this->rows, this->Ap,this->Ai, &this->Common) ;
	  }
	  
	  //change RHS to MWrap
	  BMatrix::MWrap<double>* MWrap_RHS = new  BMatrix::MWrap<double>[RHS.get_number_of_rows()];
	  for(int i=0; i<RHS.get_number_of_rows(); i++){
	      MWrap_RHS[i] = RHS.get_pointer(i,0);
	  }


	  //Do LU factorization if not done before
	  if(this->structure_has_changed){
	      this->Numeric  = klu_B_factor ( this->Ap, this->Ai, this->MWrap_Ax, this->Symbolic, &this->Common ) ;
              this->structure_has_changed = false;
	  }else if(this->values_have_changed){
	      klu_B_refactor ( this->Ap, this->Ai, this->MWrap_Ax, this->Symbolic, (klu_B_numeric*)this->Numeric, &this->Common ) ;
	      this->values_have_changed = false;
	  }
	  
	  
	  
	  //Now do the F/B substitution
	  int result = klu_B_solve(this->Symbolic, (klu_B_numeric*)this->Numeric, this->rows, Nrhs, MWrap_RHS, &this->Common);	  
	  if(!result){
		std::cerr<<"Cannot do F/B substitution";
	  }

	  sparse_free_numeric();
      }
      
      void sparse_free_numeric(){
	  klu_B_numeric* my_numeric = (klu_B_numeric*)this->Numeric ;
	  klu_B_free_numeric(&my_numeric ,  &(this->Common));
      }
      
      SBase<T>* solve(const SBase<T> &B) const{
		throw std::runtime_error("Please, code me solve(const SBase<T> &B)");
      }
      
      Sparse<T > add(const Sparse<T > &B) const{
		if(this->rows!=B.rows || this->cols!=B.cols){
			throw std::runtime_error("Can not add two sparse matrices with different sizes");
    		}
	
	  
		Sparse<T > result = *this;
		
		if(B.nnz>0){
		    typename std::list<SparseElement>::iterator B_cols_lists_iter;
		
		    for (int i = 0; i < B.cols; i++) {
			for(B_cols_lists_iter = B.cols_lists[i].begin(); B_cols_lists_iter!=B.cols_lists[i].end(); B_cols_lists_iter++){
			    result.add_to_entry(B_cols_lists_iter->row , i , B_cols_lists_iter->value);
			}
		    }
		}
	
    		
    		return result;
      }

      Sparse<T > subtract(const Sparse<T > &B) const{
		if(this->cols!=B.cols){
			throw std::runtime_error("Can not add two sparse matrices with different sizes");
    		}
	
    		Sparse<T > result = *this;
		
		if(B.nnz>0){
		    typename std::list<SparseElement>::iterator B_cols_lists_iter;
		
		    for (int i = 0; i < B.cols; i++) {
			for(B_cols_lists_iter = B.cols_lists[i].begin(); B_cols_lists_iter!=B.cols_lists[i].end(); B_cols_lists_iter++){
			    result.add_to_entry(B_cols_lists_iter->row , i , -1*(B_cols_lists_iter->value));
			}
		    }
		}
  
    		return result;
      }
      

      Sparse<T >& operator-= (const Sparse<T > &A){
	    if(this->rows!=A.rows || this->cols!=A.cols){
		throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    }
	
		
		if(A.nnz>0){
		    typename std::list<SparseElement>::iterator A_cols_lists_iter;
		
		    for (int i = 0; i < A.cols; i++) {
			for(A_cols_lists_iter = A.cols_lists[i].begin(); A_cols_lists_iter!=A.cols_lists[i].end(); A_cols_lists_iter++){
			    add_to_entry(A_cols_lists_iter->row , i , -1*(A_cols_lists_iter->value));
			}
		    }
		}
  
	    return *this;
      }

       SBase<T >& operator -=(const SBase<T > &A){
  
		const Sparse<T >* BB = dynamic_cast<const Sparse<T >*>(&A);

  		if (BB) {
			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
    			return *this-=*BB;
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
	}
      
      Sparse<T >& operator+=(const Sparse<T > &A){
	    if(this->rows!=A.rows || this->cols!=A.cols){
		throw std::runtime_error("Can not add two sparse matrices with different sizes");
	    }
	
	    if(A.nnz>0){  
		typename std::list<SparseElement>::iterator A_cols_lists_iter;

		for (int i = 0; i < A.cols; i++) {
		    for(A_cols_lists_iter=A.cols_lists[i].begin();A_cols_lists_iter!=A.cols_lists[i].end(); A_cols_lists_iter++) {
			  add_to_entry(A_cols_lists_iter->row , i, A_cols_lists_iter->value );
		    }  
		}
	    }
  
	    return *this;
      }
     
      SBase<T >& operator +=(const SBase<T > &A){
  
		const Sparse<T >* BB = dynamic_cast<const Sparse<T >*>(&A);

  		if (BB) {
			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
    			return *this+=*BB;
  		}else{
    			throw std::runtime_error("I can only handle addition of Sparse+Sparse");
  		}
	}
      
      Sparse<T > operator+ (const Sparse<T > &B) const{
		return add(B);
      }

      SBase<T >* operator+(SBase<T > const &B) const{
  		const Sparse<T >* BB = dynamic_cast<const Sparse<T >*>(&B);

  		if (BB) {
    			return new Sparse<T >( *this + *BB );
  		}else{
    			throw std::runtime_error("I can only handle addition of Sparse+Sparse");
  		}
	}

      Sparse<T > operator- (const Sparse<T > &B) const{
		return subtract(B);
      }

      SBase<T >* operator-(SBase<T > const &B) const{
  		const Sparse<T >* BB = dynamic_cast<const Sparse<T >*>(&B);

  		if (BB) {
    			return new Sparse<T >( *this + *BB  );
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
      }
	  
      Sparse<T > operator* (const Sparse<T > &B) const{
		if (this->cols != B.rows) {
      			throw std::runtime_error("Can not multibly two sparse matrices with different rows and columns");
    		}

    		int m = this->rows;
    		int n = B.cols;
 
    		Sparse<T > result(m, n);

		typename std::list<SparseElement>::iterator row_iterator , row_iterator_B;
		
    		for (int i = 0; i < B.cols; i++) {
		        for(row_iterator_B= B.cols_lists[i].begin() ; row_iterator_B!=B.cols_lists[i].end(); row_iterator_B++){
				int j = row_iterator_B->row;
				
				for (row_iterator = this->cols_lists[j].begin(); row_iterator != this->cols_lists[j].end(); row_iterator++) {
          				result.add_to_entry(row_iterator->row, i, (row_iterator->value * row_iterator_B->value));
        			}
      			}
    		}

    		return result; 
      }

     

	Dense<T > operator * (const Dense<T >& B) const{
	        if (this->cols != B.get_number_of_rows()) {
			throw std::runtime_error("Can not multible Sparse*Dense with inconsistence sizes");
		}
		
		Dense<T > result (this->rows, B.get_number_of_cols());

		typename std::list<SparseElement>::iterator row_iterator;

		
		for (int j = 0; j < B.get_number_of_cols(); j++) {	    
			    for(int col=0; col< this->cols; col++){
			        for(row_iterator= this->cols_lists[col].begin() ; row_iterator!=this->cols_lists[col].end(); row_iterator++){
					result.add_to_entry (row_iterator->row ,j, row_iterator->value * B.get(col,j) );
				}
			    }
		}
		

		return result;
	}
	
	
	Sparse<T >& operator *=(const Sparse<T > &A){
		
	
		return *this;
	}

       SBase<T >& operator *=(const SBase<T > &A){
		
  
		const Sparse<T >* BB = dynamic_cast<const Sparse<T >*>(&A);

  		if (BB) {
    			if(this->rows!=BB->rows || this->cols!=BB->cols){
				throw std::runtime_error("Can not subtract two sparse matrices with different sizes");
	    		}
			return *this *= *BB;
  		}else{
    			throw std::runtime_error("I can only handle subtraction of Sparse-Sparse");
  		}
	}
	
	Sparse<T > operator/(T val)const{
	    
	    Sparse<T > result(this->rows, this->cols);
	  
	      typename std::list<SparseElement>::iterator row_iterator;

	      for(int col=0; col< this->cols; col++){
	    	  for(row_iterator= this->cols_lists[col].begin() ; row_iterator!=this->cols_lists[col].end(); row_iterator++){
			
			SparseElement new_element;
			new_element.row = row_iterator->row;
			new_element.value = row_iterator->value/val;

			result.cols_lists[col].push_back( new_element);
			result.nnz++;
	    	  }
	      }
	    
	    return result;
	}

          Sparse<T > operator*(T val)const{
	    
	      Sparse<T > result(this->rows, this->cols);
	  
	      typename std::list<SparseElement>::iterator row_iterator;

	      if(this->nnz>0){
		  for(int col=0; col< this->cols; col++){
		      for(row_iterator= this->cols_lists[col].begin() ; row_iterator!=this->cols_lists[col].end(); row_iterator++){
			
			    SparseElement new_element;
			    new_element.row = row_iterator->row;
			    new_element.value = row_iterator->value*val;

			    result.cols_lists[col].push_back( new_element);
			    result.nnz++;
		      }
		  }
	      }

	      return result;
	  }
     
	 

	  Dense<T > operator + (BMatrix::Dense<T >& A){
		  if (A.get_number_of_rows() != this->rows
		      || A.get_number_of_cols() != this->cols) {
			    throw std::runtime_error("Can not Add Dense+Sparse with different sizes");
		  }

		  
		  BMatrix::Dense<T > result = A;
    
                  typename std::list<SparseElement>::iterator row_iterator;
		  
		    
		  for (int col = 0; col < this->cols; col++) {
		      for(row_iterator=this->cols_lists[col].begin(); row_iterator!=this->cols_lists[col].end(); col++){
			  result.add_to_entry(row_iterator->row , col, row_iterator->value );
		      }
		  }
		  	      
		  
		  return result;
	  }

	  void report_timing(){
		std::cout << "Size= "<< this->rows << "x" << this->cols<<std::endl;
		std::cout<<"Number of klu_analyze done on this matrix = "<<number_of_klu_analyze <<std::endl;
		std::cout<<"Time to do all klu_analyze= "<<time_to_do_klu_analyze << std::endl;
      	        
		std::cout<<"Number of klu_fatcor done on this matrix = "<<number_of_klu_factor <<std::endl;
		std::cout<<"Time to do all klu_factor= "<<time_to_do_klu_factor <<std::endl;

		std::cout<<"Number of klu_refatcor done on this matrix = "<<number_of_klu_refactor <<std::endl;
		std::cout<<"Time to do all klu_refactor= "<<time_to_do_klu_refactor <<std::endl;
	  }

	 klu_symbolic* get_Symbolic(){return this->Symbolic;}
	 klu_numeric* get_Numeric(){return this->Numeric;}

	 void use_Symbolic(klu_symbolic* _symbolic){
		this->Symbolic = _symbolic;
		this->structure_has_changed = false;
	 }
	 void use_Numeric(klu_numeric* _numeric){
		this->Numeric = _numeric;
		this->structure_has_changed = false;
	 }

      protected:
      void runtime_error(const char* arg1);

      template<class U>
	friend std::ostream& operator << (std::ostream &out ,  Sparse<U >& B);
};
  
  

  
 
template<class U> 
std::ostream& operator << (std::ostream &out ,  Sparse<U>& B){

      /*const_cast<Sparse<U>&>(B).create_ccs();
      
      out<<"this->this->Ap=[";
      for(int i=0; i<B.cols+1; i++){
	  out<<B.Ap[i]<<" ";
      }
      out<<"]"<<std::endl;
      
      out<<"this->this->Ai=[";
      for(int i=0; i<B.nnz; i++){
	  out<<B.Ai[i]<<" ";
      }
      out<<"]"<<std::endl;
      
      out<<"this->this->Ax=[";
      for(int i=0; i<B.nnz; i++){
	  out<<(B.Ax[i])<<" ";
      }
      out<<"]"<<std::endl;*/
      
      for(int i=0; i<B.cols; i++ ){
	  for(auto e : B.cols_lists[i]){
	      out << e.row<<" "<<i<<" "<<e.value <<std::endl;
	  }
      }
      
      return out;
} 
  
};

#endif