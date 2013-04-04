//
// C++ Implementation: subcircuit
//
// Description: 
//
//
// Author: Mina Farhan <madel@diode.doe.carleton.ca>, (C) 2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "subcircuit.h"
#include <string>



/*-------------------SubCircuitInstance Class functions------------------------*/

///Create an instance of this subicrcuit with the name "name"
SubCircuitInstance* SubCircuit::create_instance(std::string name , std::vector<std::string> terminals){
      return new SubCircuitInstance(name, terminals, this);
}


void SubCircuitInstance::add_my_nodes(Circuit* circuit , const std::vector<std::string>& append_to_node_name){
     std::vector<Element*>::iterator iter;     
     
      
     
     //check if there is any component that needs to add an extra variable
     for(iter=parent_subcircuit->components.begin(); iter!=parent_subcircuit->components.end(); iter++){
	  
         //the vector of strings to append to node names
	  std::vector<std::string> extra_string;
     
	  //We have to check if this node is from the terminals...If it is a terminal node then we do not append the name of the subcircuit to it
	  std::vector<std::string> element_terminals = (*iter)->get_terminals_name();
	  for(int i=0; i<element_terminals.size(); i++){
	      //if this element has a node that is a terminal, then do not append something to it
	      if(find(parent_subcircuit->terminals.begin() , parent_subcircuit->terminals.end() , element_terminals[i]) !=  parent_subcircuit->terminals.end() ){
		  extra_string.push_back(append_to_node_name + name + ".");
	      }else{
		  extra_string.push_back("");
	      }
	  }
	  
	  //Add the nodes of this element
	  (*iter)->add_my_nodes(circuit , extra_string); 
     } 
  
}

void SubCircuitInstance::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
      
}
