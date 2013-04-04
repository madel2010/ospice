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
	std::vector<std::string> terminals;
	std::string name;

public:
	SubCircuit(std::string _name, std::vector<std::string> _terminals):name(_name), terminals(_terminals){};
	
	///Create an instance of this subicrcuit with the name "name"
	inline SubCircuitInstance* create_instance(std::string name , std::vector<std::string> terminals);
	  
	
};


class SubCircuitInstance: public Element{
  
private:
    SubCircuit* parent_subcircuit;
    std::vector<std::string> terminals;
    
    int num_nodes_to_add; //This is the number of nodes that the subcircuit needs. 
			 //Note: this number should execlude the terminals nodes as they are already in the main circuit
    
public:
    SubCircuitInstance() {}
    
    SubCircuitInstance(std::string _name,  std::vector<std::string> _terminals, SubCircuit* _subcircuit):terminals(_terminals){
	name = _name;
	parent_subcircuit = _subcircuit;
	
	num_nodes_to_add = 0;
    }
  
    /////Start:Must be defined as it it inherited from Element Base class
    void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
    inline bool is_linear(){return parent_subcircuit->is_linear;}
    void add_my_nodes(Circuit* circuit , const std::vector<std::string>& append_to_node_name); 
    /////END:Must be defined as it it inherited from Element Base class

    //Returns the names of the terminals 
    std::vector<std::string> get_terminals_name(){    
	    return terminals;
    }
};

#endif
