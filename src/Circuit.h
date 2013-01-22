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


#ifndef CIRCUIT_H
#define CIRCUIT_H


#include "BMatrix.h"
#include "element.h"
#include <vector>
#include <map>



class Element;
class NonLinElement;
class Probe;
class Source;
class Analysis;
class DC;

class Circuit
{
  
private:
     BMatrix::Sparse<double> G;
     BMatrix::Sparse<double> C;
     BMatrix::Sparse<double> fx; //the non linear vector
     BMatrix::Dense<double> B; 
     
     
     //Vector to hold list of elemets
     std::vector<Element*> components;
     
     //Vector to hold list of elemets
     std::vector<Source*> sources;
     
     //Vector to hold list of non linear elements. This is used to update the vector fx each iteration
     std::vector<NonLinElement*> Non_Linear_Elements;
     
     //Vector to hold list of analysis to be done on the circuit
     std::vector<Analysis*> required_analysis;
     
     //associative array that maps the node name/current_name to its index
     std::map< std::string , int> mna_variable_indices;
     
     //Vector to hold list of elemets
     std::vector<Probe*> Probes;
     
     //associative array that saves the elements that addes extra nodes for currents like inductors. 
     //This is usefeull when we need to add mutual inductance. 
     //In this case we do not need to search in all elements
     // map<name , index_of_current>
     std::map< std::string , int> current_elements;
     
     void attach_elements();
     
     bool is_linear;
     
     bool have_dc_solution;
     DC* DC_Analysis; //saves a pointer to the DC analysis
     
public:
     Circuit();
    ~Circuit();
    
    int get_variable_index(std::string node_name); //return the index of the node
    
    void add_mna_variable(std::string node_name);
    
    inline void add_probe(Probe* _Probe){
	  Probes.push_back(_Probe);
    }
    
    void add_current_element(std::string element_name);
    
    inline void add_source(Source* _src){
      sources.push_back(_src);
    }
    
    
    inline void add_NonLinElement(NonLinElement* E){
      Non_Linear_Elements.push_back(E);
    }

    
    void operator << (Element* E);
    
    inline void operator << (Analysis* an){
	  required_analysis.push_back(an);
    }
    
    int size_of_mna(){return mna_variable_indices.size(); }
    
    bool is_circ_linear() {return is_linear; }
    
    void update_probes(double time, const double* solution);
    
    void update_sources(double time);
    
    void update_fx(const double* solution, double time);
 
    void start_analysis();
    
    void plot_probes();
    
    const BMatrix::Dense<double>& get_dc_solution();
};

#endif // CIRCUIT_H
