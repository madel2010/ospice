/*
    Copyright (c) 2012, <copyright holder> <email>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY <copyright holder> <email> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <copyright holder> <email> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "resistor.h"
#include "Sparse.h"
#include "MatrixBase.h"

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include "ShuntingYard.h"

using boost::lexical_cast;
using boost::bad_lexical_cast;
    
/*-------------Linear Resisor--------------*/
void resistor::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    
    if(n1_index > -1){
	G.add_to_entry(n1_index , n1_index , 1/value);
    }
    
    if(n2_index > -1) {
        G.add_to_entry(n2_index , n2_index , 1/value);
    }
    
    if(n2_index > -1 && n1_index >-1 ) {
	G.add_to_entry(n1_index , n2_index , -1/value);
	G.add_to_entry(n2_index , n1_index , -1/value);
    }
}


void resistor::add_my_nodes(Circuit* circuit){
    n1_index = circuit->add_mna_variable(n1);
    n2_index = circuit->add_mna_variable(n2);
}



/*-------------Non Linear Resisor--------------*/
nonlin_resistor::nonlin_resistor(std::string _n1, std::string _n2, const char* _curr_expression):TwoTerminal(_n1,_n2){

	     name = std::string("nlR")+".+"+n1+".-"+n2;

             //Remove all spaces from the expression to make it more compact when parsing it	     
	     Expression =  boost::regex_replace(std::string(_curr_expression),boost::regex("\\s+"),"");
	     
	     //add paranthesis to the expression. It is required for the  shunting_yard algorithm to make sure the expression is done.

             Expression = std::string("(") + Expression;
	     Expression += ")";
	     
	}
	
nonlin_resistor::nonlin_resistor(std::string _name, std::string _n1, std::string _n2, const char* _curr_expression):TwoTerminal(_n1,_n2){
	
	     name =_name;
	     
	     Expression = boost::regex_replace(std::string(_curr_expression),boost::regex("\\s+"),"");
	     
	     //add paranthesis to the expression. It is required for the  shunting_yard algorithm to make sure the expression is done
	     Expression = std::string("(") + Expression;
	     Expression += ")";
	     
	};

void nonlin_resistor::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    
    circ->add_NonLinElement(this);
    
    //call the shunt_yard algorithm to parse the nonlinear expression
    F = shunting_yard(Expression , expression_depends_on , circ);
    
    //Find dF/dx. We assume that all the tree nodes has been evaluated in the shunting_yard algorithm, and F has been evaluated
    std::vector< std::pair<Symbolic , int> >::iterator iter;
    for(iter=expression_depends_on.begin(); iter!=expression_depends_on.end(); iter++){
	dFdx.push_back(df(F,iter->first)); //deffierentaite F with respect to each variable it depends on in expression_depends_on
    }
}


void nonlin_resistor::add_my_nodes(Circuit* circuit ){
    n1_index = circuit->add_mna_variable(n1);
    n2_index = circuit->add_mna_variable(n2);
}




void nonlin_resistor::update_fx(BMatrix::Dense<double>& fx, const double* solution){
      //firxt Evalue the expression at the specific time
      
      
      std::vector< std::pair<Symbolic , int> >::iterator iter;
      double value;
      Symbolic temp_F = F;

      for(iter=expression_depends_on.begin(); iter!=expression_depends_on.end(); iter++){
	  //get the value of the voltages that the expression depends on
	  value = solution[iter->second];
	  
	  temp_F = temp_F.subst(iter->first, value);
      }



      try
      {
	  value = lexical_cast<double>(temp_F);
      }
      catch(boost::bad_lexical_cast& e)
      {
	  throw std::runtime_error("Expression can not be converted to a value, something wrong");
      }
    
      if(n1_index > -1) fx.put(n1_index,0,value);
      if(n2_index > -1) fx.put(n2_index,0,-value);
}

//This function gives the result of J by updating the entries of the jacobian matrix 
void nonlin_resistor::update_J(BMatrix::Sparse<double>& J, const double* solution){
  
      std::vector< std::pair<Symbolic , int> >::iterator iter_depends_on;
      std::vector<Symbolic>::iterator iter_dFdx;
      
      double value;
      Symbolic temp;
      int dfdx_counter=0;
      
      for( iter_dFdx=dFdx.begin();  iter_dFdx!=dFdx.end(); iter_dFdx++){
	  
	  
	  temp = (*iter_dFdx);

	  for(iter_depends_on=expression_depends_on.begin(); iter_depends_on!=expression_depends_on.end(); iter_depends_on++){
	      temp = temp.subst(iter_depends_on->first , solution[iter_depends_on->second]);
	  }

	  try
	  {
	      value = lexical_cast<double>(temp);
	  }
	  catch(boost::bad_lexical_cast& e)
	  {
	      throw std::runtime_error("DF/DX can not be converted to a value, something wrong");
	  }
	  
	  if(n1_index > -1){
	      J.put(n1_index , expression_depends_on[dfdx_counter].second , value);
	  }
    
	  if(n2_index > -1) {
	      J.put(n2_index , expression_depends_on[dfdx_counter].second , value);
	  }
      
	  dfdx_counter++;
      }
      
}
