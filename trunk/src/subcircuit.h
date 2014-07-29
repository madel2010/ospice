//
// C++ Interface: subcircuit
//
// Description: 
//
//
// Author: Mina Farhan <madel@diode.doe.carleton.ca>, (C) 2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SUBCIRCUIT_H
#define SUBCIRCUIT_H

#include <vector>
#include <string>
#include "Circuit.h"
#include "element.h"

class SubCircuitInstance;

class SubCircuit : public Circuit
{

friend class SubCircuitInstance;
  
private:
	std::string name;
        std::vector<std::string> terminals;
	
public:
	SubCircuit(std::string _name, std::vector<std::string> _terminals):name(_name), terminals(_terminals){};
	
	///Create an instance of this subicrcuit with the name "name"
	SubCircuitInstance* create_instance(std::string name , std::vector<std::string>& terminals);
	
	std::string get_name() const{return name;}
	  
	
};


class SubCircuitInstance: public Element{
  
private:
    SubCircuit* my_subcircuit;
    std::vector<std::string> terminals;
    
    std::string prepend_to_nodes; //a string that would be prepended to all nodes of this instance
    
    std::list<Element*> cloned_elements;
    
public:
    SubCircuitInstance() {}
    ~SubCircuitInstance();
    
    SubCircuitInstance(std::string _name,  std::vector<std::string>& _terminals, SubCircuit* _subcircuit):terminals(_terminals){
	name = _name;
	my_subcircuit = _subcircuit;
	prepend_to_nodes = "";
    }
  
    SubCircuitInstance* clone(){ return new SubCircuitInstance(*this); }
     
    /////Start:Must be defined as it it inherited from Element Base class
    void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
    
    inline bool is_linear(){return my_subcircuit->is_linear;}
    
    //The Subcircuit does not add currents to the MNA, then return -1
    int is_current_element(){return -1;}
	
    void add_my_nodes(Circuit* circuit); 
    /////END:Must be defined as it it inherited from Element Base class

    //Returns the names of the terminals 
    std::vector<std::string> get_terminals(){    
	    return terminals;
    }
    
    void prepend_name(std::string p){name = p + name;}
    
    void prepend_nodes(std::string p){ prepend_to_nodes = p; }
};

#endif
