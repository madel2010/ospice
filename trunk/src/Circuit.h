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
#include "symbolicc++.h"


class Element;
class NonLinElement;
class Probe;
class Source;
class Analysis;
class DC;



class Circuit
{


protected:
     bool is_linear;

     
     //Vector to hold list of elemets
     std::list<Element*> components;
     
     //Vector to hold list of elemets
     std::list<Source*> sources;
     
     //Vector to hold list of non linear elements. This is used to update the vector fx each iteration
     std::list<NonLinElement*> Non_Linear_Elements;
     
     //Map for Parameters specified in the netlist map[Name] = expression;
     std::map<std::string , std::string> Parameters;

     //associative array that saves the elements that addes extra nodes for currents like inductors. 
     //This is usefeull when we need to add mutual inductance. 
     //In this case we do not need to search in all elements
     // map<name , <index_of_current,value of inductor> >
     std::map< std::string , std::pair<int,double> > inductors;


private:
     BMatrix::Sparse<double> G;
     BMatrix::Sparse<double> C;
     BMatrix::Sparse<double> J;
     BMatrix::Dense<double> fx; //the non linear vector
     BMatrix::Dense<double> B; 
     
     //associative array that maps the node name/current_name to its index
     std::map< std::string , int> mna_variable_indices;
     
     //associative array that saves if there are nodes have different names but same inex. very usefull in subcircuits
     std::map< std::string , std::string> similar_nodes;
     
     //Vector to hold list of elemets
     std::list<Probe*> Probes;

     //Vector to hold list of analysis to be done on the circuit
     std::vector<Analysis*> required_analysis;
     
     bool have_dc_solution;
     DC* DC_Analysis; //saves a pointer to the DC analysis
 
     virtual void attach_elements();
public:
     Circuit();
    ~Circuit();
    
    //returns node_index if exists, -1 if ground , -2 if not exists
    int get_variable_index(std::string var_name); //return the index of the node
    
    int add_mna_variable(std::string var_name);
    
    //This makes two node have the same index.
    void alias_two_nodes(std::string node1 , std::string node2 ){similar_nodes[node1] = node2;}
    
    
    void add_inductor_index(std::string inductor_name, double value);
    
    //this function gets the index of the current of the current_element
    std::pair<int,double> get_inductor_element_index(std::string element_name){
	return inductors[element_name];
    }
    
    inline void add_source(Source* _src){
      sources.push_back(_src);
    }

    void add_parameter(std::string name, std::string expression);
    bool get_parameter_expression(std::string name, std::string& return_expression); //return false if parametere does not exist

    inline void add_NonLinElement(NonLinElement* E){
      Non_Linear_Elements.push_back(E);
    }
    
     inline void add_probe(Probe* _Probe){
	  Probes.push_back(_Probe);
     }

     void operator << (Element* E);

    inline void operator << (Analysis* an){
	  required_analysis.push_back(an);
    }
    
    Element* search_elements(std::string name);
    
    int size_of_mna(){return mna_variable_indices.size(); }
    
    void update_probes(double time, const double* solution);
    
    void update_sources(double time);
    
    void update_fx( const double* solution);
    
    void update_J(const double* solution);
 
    void start_analysis();
    
    void plot_probes();
    
    const BMatrix::Dense<double>& get_dc_solution();

    bool is_circ_linear() {return is_linear; }
};

#endif // CIRCUIT_H
