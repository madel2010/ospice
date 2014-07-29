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


#ifndef SOURCE_H
#define SOURCE_H

#include <element.h>
#include <string>
#include "SourceFunc.h"
#include "Circuit.h"

class Source : public TwoTerminal
{

protected:
  SourceFunc* Func;
  
public:
    Source(std::string _n1, std::string _n2):TwoTerminal(_n1,_n2){};
    ~Source(){
	delete Func;
    };
    virtual void update_B(BMatrix::Dense<double> &B,double time)=0;
};


/*-------CurrentSource Class ----------------- */
class CurrentSource : public Source
{
private:
      
public:
     CurrentSource(std::string _n1, std::string _n2, SourceFunc* _Func):Source(_n1,_n2){
	  name = std::string("CurrentSource")+".+"+n1+".-"+n2;
	  Func = _Func;
     };
     
     CurrentSource(std::string _name, std::string _n1, std::string _n2, SourceFunc* _Func):Source(_n1,_n2){
	  Func = _Func;
	  name = _name;
     };
     
    
    CurrentSource* clone(){ return new CurrentSource(*this); }
     
    void update_B(BMatrix::Dense<double> &B, double time);
    
    ///add the requird nodes to the main circuit
    void add_my_nodes(Circuit* circuit);
    
    //Returns the names of the terminals 
    std::vector<std::string> get_terminals_name(){
	    std::vector<std::string> result;
	    result.push_back(n1);
	    result.push_back(n2);
	    
	    return result;
    }
	
    void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
    
    bool is_linear(){ return true;} ;
    
    //The current source does not add currents to the MNA, then return -1
    int is_current_element(){return -1;}
    
};

/*-------Voltage Source Class ----------------- */
class VoltageSource : public Source
{
private:

      int current_index; //the index of the extra current variable
      
      
public:
     VoltageSource(std::string _n1, std::string _n2, SourceFunc* _Func):Source(_n1,_n2){
	  name = std::string("VoltageSource")+".+"+n1+".-"+n2;
	  Func = _Func;
     };
     
     VoltageSource(std::string _name, std::string _n1, std::string _n2, SourceFunc* _Func):Source(_n1,_n2){
	  Func = _Func;
	  name = _name;
     };
     
    VoltageSource* clone(){ return new VoltageSource(*this); }
    
    void update_B(BMatrix::Dense<double> &B, double time);
    
    ///add the requird nodes to the main circuit
    void add_my_nodes(Circuit* circuit);
    
    void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
    
    //Returns the names of the terminals 
    std::vector<std::string> get_terminals_name(){
	    std::vector<std::string> result;
	    result.push_back(n1);
	    result.push_back(n2);
	    
	    //the current
	    result.push_back(name + ".I");
	    
	    return result;
    }
    
    bool is_linear(){ return true;} ;
    
    //The Voltage source adds currents to the MNA, then return the index of its current index
    int is_current_element(){return current_index;}
    
};

#endif // SOURCE_H
