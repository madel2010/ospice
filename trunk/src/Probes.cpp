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

#include "Circuit.h"
#include "Probes.h"
#include "plot.h"

/*-----------------Voltage porbles------------*/
VoltageProbe::VoltageProbe(std::string _name , std::string _n1, std::string _n2){
    n1 = _n1;
    n2 = _n2;
    name = _name;
}


VoltageProbe::VoltageProbe(TwoTerminal* element){
    n1 = element->get_n1();
    n2 = element->get_n2();
    name = element->get_name();
  
}

void VoltageProbe::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C , Circuit* circ){
    n1_index = circ->get_variable_index(n1);
    n2_index = circ->get_variable_index(n2);
    
    circ->add_probe(this);
    
}

///This function gets the solution at each time step from the Circuit class and put it in the probes vectors
void VoltageProbe::get_data(double time, const double* solution){
      double value;
      
      if(n1_index >-1){ //>-1 means it is not ground, becasue gnd_index = -1
	  value = solution[n1_index];
      }
      
      if(n2_index >-1){ //>-1 means it is not ground, becasue gnd_index = -1
	  value -= solution[n2_index];
      }
      
      data.push_back(value);
      time_points.push_back(time);
}
 
///THis function plots the porbe data using gnuplot
void VoltageProbe::plot(){
    if(data.size() != time_points.size()){
	throw std::runtime_error("Something wrong happened. A probe doesnot have the same number of data and time points");
    }
    
    if(data.size() == 0 && time_points.size() == 0){
	std::string err = std::string("Probe with name '") + name + std::string("' Does not have any data"); 
	std::cout<< err;
    }
    
    
    my_plot.plot( time_points, data, name); //this function is defined in plot.h to plot the data of the probe
}