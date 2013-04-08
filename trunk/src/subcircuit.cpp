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
#include <algorithm> 


/*-------------------SubCircuitInstance Class functions------------------------*/

///Create an instance of this subicrcuit with the name "name"
SubCircuitInstance* SubCircuit::create_instance(std::string name , std::vector<std::string> terminals){
      return new SubCircuitInstance(name, terminals, this);
}


void SubCircuitInstance::add_my_nodes(Circuit* circuit){
     
      //Alias the terminal nodes to be the same as the nodes it is attached to
      for(int i=0; i<my_subcircuit->terminals.size(); i++){
	  std::string node_name = name + my_subcircuit->terminals[i];
	  circuit->alias_two_nodes( node_name , this->terminals[i]);
      }
      
     //Add the elements to the main circuit with new names and nodes reflectin the subcircuit instance
     std::list<Element*>::iterator iter;     
     for(iter=my_subcircuit->components.begin(); iter!=my_subcircuit->components.end(); iter++){
         //just clone the elements to a new element with same value but different node names and attach it to the circuit using <<
         Element* newelem = (*iter)->clone();
	 newelem->prepend_name(name+".");
	 newelem->prepend_nodes(name+".");
         cloned_elements.push_back(newelem);
     } 
  
}

void SubCircuitInstance::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
      std::list<Element*>::iterator iter; 
      for(iter=cloned_elements.begin(); iter!=cloned_elements.end(); iter++){
	  (*iter)->write_stamp(G,C,circ);  
      }
}

SubCircuitInstance::~SubCircuitInstance(){
     //do not delete them because they will be deleted from circuit class
     //std::list<Element*>::iterator iter;     
     //for(iter=cloned_elements.begin(); iter!=cloned_elements.end(); iter++){
     //    delete (*iter);
     //} 
}
