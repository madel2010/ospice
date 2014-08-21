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


#ifndef ELEMENT_H
#define ELEMENT_H

#include "Sparse.h"
#include "MatrixBase.h"
#include "Circuit.h"
#include <boost/regex.hpp>

class Circuit;

class Element
{
  
protected:
      std::string name; //The name of the element
      
      //this is if we need to sort the elements of the circuit. 
      //For example, the probes should be sorted at the end, so we give it index 2
      //Any element by default = 1; probes = 2
      int element_order_index; 
      
public:
    //base constructor is called first
    Element(){element_order_index=1;};
    
    virtual Element* clone()=0;
    
    virtual void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ)=0;
    virtual bool is_linear()=0;
    virtual int is_current_element()=0; //does this element add a current to the MNA. 
					//If No, it should return -1, Otherwise it should return the index of the current element
    
    ///add the node name to the circuit, 
    virtual void add_my_nodes(Circuit* circuit)=0; 
    
    ///get the name of the element
    virtual std::string get_name(){return name;}
    
    //get the order of this element in the circuit_elements list
    int get_order_index(){return element_order_index;}
    
    virtual void prepend_name(std::string p){name = p + name;}
    
    virtual void prepend_nodes(std::string p)=0;
    
    ///Get the terminals of this element. 
    ///It returns the names of the terminals of this element
    virtual std::vector<std::string> get_terminals()=0;
    
    virtual ~Element(){}
};

class NonLinElement
{
protected:
    std::string Expression;
    
public:
    NonLinElement(){};
    virtual void update_fx(BMatrix::Dense<double>& fx, const double* solution)=0;
    virtual void update_J(BMatrix::Sparse<double>& J, const double* solution)=0;
    
    //add p before the node or currents in thenonlinear expression. Used when we have subcircuits
    virtual void prepend_expression(std::string p){
      boost::regex e ("([vV])\\([ \t]*([A-Za-z_]?[0-9]*)[ \t]*\\)");   
      Expression = boost::regex_replace (Expression,e,std::string("$1(")+p+"$2)");
    }
    
    virtual ~NonLinElement(){}
    
};

/*------------------------------------------*/
///Two Terminal Element Class
class TwoTerminal : public Element
{
friend class CurrentProbe;

protected:
    std::string n1; //First node name
    std::string n2; //second node name
    
    int n1_index; //the index of the first node
    int n2_index; //the index of the second node
    
public:
  
    TwoTerminal(){};
    TwoTerminal(std::string _n1, std::string _n2){
	n1 = _n1;
	n2 = _n2;
    }
    
    std::string get_n1(){return n1;}
    std::string get_n2(){return n2;}

    virtual ~TwoTerminal(){}

    virtual std::vector<std::string> get_terminals(){
	std::vector <std::string> result;
	result.push_back(n1);
	result.push_back(n2);
	return result;
    }
    
    //add p before the node names. Used when we have subcircuits
    virtual void prepend_nodes(std::string p){
	n1 = p + n1;
	n2 = p + n2;
    }
};

/*------------------------------------------*/
///Four Terminal Element Class
class FourTerminal : public Element
{

protected:
    std::string in1; //First input node name
    std::string in2; //second input node name
    
    std::string out1; //First output node name
    std::string out2; //second output node name
    
    int in1_index; //the index of the first input node
    int in2_index; //the index of the second input node
    
    int out1_index; //the index of the first output node
    int out2_index; //the index of the second output node
    
public:
  
    FourTerminal(){};
    FourTerminal(std::string _in1, std::string _in2, std::string _out1, std::string _out2 ){
	in1 = _in1;
	in2 = _in2;
	out1 = _out1;
	out2 = _out2;
    }
    
    std::string get_in1(){return in1;}
    std::string get_in2(){return in2;}
    std::string get_out1(){return out1;}
    std::string get_out2(){return out2;}

    virtual ~FourTerminal(){}

    virtual std::vector<std::string> get_terminals(){
	std::vector <std::string> result{in1 , in2, out1, out2};
	return result;
    }
    
    virtual void prepend_nodes(std::string p){
	in1 = p + in1;
	in2 = p + in2;
	out1 = p + out1;
	out2 = p + out2;
    }
};

#endif // ELEMENT_H
