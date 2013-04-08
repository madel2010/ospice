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


#ifndef INDUCTOR_H
#define INDUCTOR_H

#include "element.h"


class Inductor : public TwoTerminal
{
private:
	double value;
	int current_index;
	
public:
		
	Inductor(std::string _n1, std::string _n2, double _value):TwoTerminal(_n1,_n2),value(_value){
	     name = std::string("L")+".+"+n1+".-"+n2;
	}
	
	Inductor(std::string _name, std::string _n1, std::string _n2, double _value):TwoTerminal(_n1,_n2),value(_value){
	  name = _name;
	};
	
	Inductor* clone(){return new Inductor(*this);} 
	
	bool is_linear(){return true;}
	
	virtual void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
	
	///add the requird nodes to the main circuit
	void add_my_nodes(Circuit* circuit);
	
	std::vector<std::string> get_terminals_name(){
	    std::vector<std::string> result;
	    result.push_back(n1);
	    result.push_back(n2);
	    
	    //the current_index
	    result.push_back(name + ".I");
	    
	    return result;
	}
	
};


class MututalInductor : public TwoTerminal
{
private:
	double value;

	
public:
		
	MututalInductor(std::string _n1, std::string _n2, double _value):TwoTerminal(_n1,_n2),value(_value){
	     name = std::string("K")+".+"+n1+".-"+n2;
	}
	
	MututalInductor(std::string _name, std::string _n1, std::string _n2, double _value):TwoTerminal(_n1,_n2),value(_value){
	  name = _name;
	};
	
	
	bool is_linear(){return true;}
	
	virtual void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
	
	///add the requird nodes to the main circuit,
	void add_my_nodes(Circuit* circuit){
	  //in case this mutual inductor is in subcircuit, then we have to change the name of the inductors to reflect the subcircuit name
	  n1  = n1;
	  n2  = n2;
	};
	
	std::vector<std::string> get_terminals_name(){
	    std::vector<std::string> result;
	    
	    return result;
	}
	
};

#endif // INDUCTOR_H
