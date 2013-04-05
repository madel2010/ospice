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
     std::vector<Element*>::iterator iter;     
     
     
     
     //check if there is any component that needs to add an extra variable
     for(iter=my_subcircuit->components.begin(); iter!=my_subcircuit->components.end(); iter++){
         //just clone the elements to a new element with same value but different node names and attach it to the circuit using <<
     } 
  
}

void SubCircuitInstance::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
      std::vector<Element*>::iterator iter; 
      for(iter=my_subcircuit->components.begin(); iter!=my_subcircuit->components.end(); iter++){
	  (*iter)->write_stamp(G,C,circ);  
      }
}
