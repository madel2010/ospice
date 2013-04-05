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

class Circuit;

class Element
{
  
protected:
      std::string name; //The name of the element

public:
    Element(){};
    virtual void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ)=0;
    virtual bool is_linear()=0;
    
    ///add the node name to the circuit, 
    virtual void add_my_nodes(Circuit* circuit)=0; 
    
    ///get the name of the element
    virtual std::string get_name(){return name;}
    
    ///Get the terminals of this element. 
    ///It returns the names of the terminals of this element
    virtual std::vector<std::string> get_terminals_name()=0;
    
    virtual ~Element(){}
};

class NonLinElement
{

public:
    NonLinElement(){};
    virtual void update_fx(BMatrix::Dense<double>& fx, const double* solution)=0;
    virtual void update_J(BMatrix::Sparse<double>& J, const double* solution)=0;
    
    virtual ~NonLinElement(){}
    
};

/*------------------------------------------*/
///Two Terminal Element Class
class TwoTerminal : public Element
{

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


};


#endif // ELEMENT_H
