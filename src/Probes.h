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


#ifndef PROBES_H
#define PROBES_H

#include "element.h"
#include "Inductor.h"
#include <vector>
#include <string>
#include "plot.h"


class Circuit;

class Probe : public TwoTerminal
{
  
protected:

    std::vector<double> data;
    std::vector<double> time_points;
    
    Plot* my_plot; //this keeps the class of the plot window of gnuplot
public:
    
    //element_order_index is defined in Element class. It gives an order to the elements to write in the MNA first. 
    //we need probes to be the last elements to make sure that all nodes/elements are added first.
    Probe(){element_order_index=3; my_plot = nullptr;};
    Probe(Plot* _my_plot){element_order_index=3; my_plot = _my_plot;};
    virtual ~Probe(){
      if(my_plot){ 
	delete my_plot;
	my_plot = nullptr;
      }
    };

    virtual void get_data(double time , const double* solution)=0;
    virtual void plot()=0;
    
};

/*---------------Voltage Probe Class -----------*/
class VoltageProbe : public Probe
{
private:
    int n1_index;
    int n2_index;
    
    
public:
    VoltageProbe(std::string _name, std::string _n1, std::string _n2, Plot* _my_plot = nullptr);
    VoltageProbe(TwoTerminal* element, Plot* _my_plot = nullptr);  
   
    VoltageProbe* clone(){ return new VoltageProbe(*this); }
    
    void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
    
    bool is_linear(){return true;}
    
    //The VoltageProbe does not add currents to the MNA, then return -1
    int is_current_element(){return -1;}
	
    ///add the requird nodes to the main circuit
    void add_my_nodes(Circuit* circuit);
    
    void get_data(double time , const double* solution);
    
    void plot();
    
    //Returns the names of the terminals 
    std::vector<std::string> get_terminals_name(){
	    std::vector<std::string> result;
	    result.push_back(n1);
	    result.push_back(n2);
	    
	    return result;
    }
};

/*---------------Current Probe Class -----------*/
class CurrentProbe : public Probe
{
  
private:
   int current_index;
   TwoTerminal* my_element;
   std::string my_element_name;
   
   ShortCircuit* SC;
   
public:
  
    CurrentProbe(std::string _name, TwoTerminal* element, Plot* _my_plot = nullptr):Probe(_my_plot){
      name = _name;
      my_element = element;
      SC = nullptr;
    }
    CurrentProbe(TwoTerminal* element, Plot* _my_plot = nullptr):Probe(_my_plot){
	name = std::string("I(")+element->name+std::string(")"); 
	SC = nullptr; 
	my_element = element;
    } 
   
    CurrentProbe(std::string _name, std::string element_name, Plot* _my_plot = nullptr):Probe(_my_plot){
      name = _name;
      my_element = nullptr;
      SC = nullptr;
      my_element_name = element_name;
    }
    
    CurrentProbe* clone(){ return new CurrentProbe(*this); }
    
    void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
    
    bool is_linear(){return true;}
    
    //The current Probe adds short circuit which requires currents added to the MNA, then return current index
    int is_current_element(){return current_index;}
	
    ///add the requird nodes to the main circuit
    void add_my_nodes(Circuit* circuit);
    
    void get_data(double time , const double* solution);
    
    void plot();
    
    
};

#endif // PROBES_H
