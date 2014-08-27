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


#ifndef RESISTOR_H
#define RESISTOR_H

#include "element.h"
#include "BMatrix.h"
#include <string>
#include "Circuit.h"
#include <map>


#include "symbolicc++.h"

class resistor : public TwoTerminal
{
  
private:
	double value;
	
	
public:
		
	resistor(std::string _n1, std::string _n2, double _value):TwoTerminal(_n1,_n2),value(_value){
	     name = std::string("R")+".+"+n1+".-"+n2;
	}
	
	resistor(std::string _name, std::string _n1, std::string _n2, double _value):TwoTerminal(_n1,_n2),value(_value){
	  name =_name;
	};

	resistor* clone(){ return new resistor(*this); }
	
	bool is_linear(){return true;}
	
	//The Resistor does not add currents to the MNA, then return -1
	int is_current_element(){return -1;}
	
	virtual void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
	
	///add the requird nodes to the main circuit
	void add_my_nodes(Circuit* circuit);
	
	//Returns the names of the terminals 
	std::vector<std::string> get_terminals_name(){
	    std::vector<std::string> result;
	    result.push_back(n1);
	    result.push_back(n2);
	    
	    return result;
	}
	
};

/*----------Class of non Linear Resistors---------------*/
class nonlin_resistor: public TwoTerminal, public NonLinElement
{
  
private:

      Symbolic F; //The non linear function 
      std::vector<Symbolic> dFdx; //The derivative of df/dx with respect to each node that F depends on
      
      std::vector< std::pair<Symbolic , int> > expression_depends_on; //The node index and references that the expression depends on. 
      

public:
		
	nonlin_resistor(std::string _n1, std::string _n2, const char* _curr_expression);
	
	nonlin_resistor(std::string _name, std::string _n1, std::string _n2, const char* _curr_expression);

	nonlin_resistor* clone(){ return new nonlin_resistor(*this); }
	
	bool is_linear(){return false;}
	
	//The Resistor does not add currents to the MNA, then return -1
	int is_current_element(){return -1;}
	
	virtual void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
	
	///add the requird nodes to the main circuit
	void add_my_nodes(Circuit* circuit);
	
	//This function adds the falue to the nonlinear vector fx of the circuit
       void update_fx(BMatrix::Dense<double>& fx, const double* solution);
      
       void update_J(BMatrix::Sparse<double>& J, const double* solution);
       
       //Returns the names of the terminals 
       std::vector<std::string> get_terminals_name(){
	    std::vector<std::string> result;
	    result.push_back(n1);
	    result.push_back(n2);
	    
	    return result;
       }
       
       //nonlinear resistor will have its own copy of prepend_nodes
       void prepend_nodes(std::string p){
	  TwoTerminal::prepend_nodes(p);
	  
	  prepend_expression(p);
       }
};

#endif // RESISTOR_H
