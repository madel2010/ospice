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
#include "Source.h"
#include "Analysis.h"

Circuit::Circuit()
{
    //thie circuit is linear until we start adding elements using <<
    is_linear = true;
}

Circuit::~Circuit()
{ 
     std::list<Element*>::iterator it_el;
     for(it_el=components.begin(); it_el!=components.end();it_el++){
	  delete (*it_el);
     } 
	    
     /* we do not need that because sources are also added in the components vectro
     std::vector<Source*>::iterator it_src;
     for(it_src=sources.begin(); it_src!=sources.end(); it_src++){
	  delete (*it_src);
     } 
     */

     std::vector<Analysis*>::iterator it_an;
     for(it_an=required_analysis.begin(); it_an!=required_analysis.end();it_an++){
	delete (*it_an);
     }     
}

int Circuit::add_mna_variable(std::string var_name){
    
    int index = get_variable_index(var_name);
    
    if(index==-1 ){ //means it is ground
	return -1;
    }else if( index == -2 ){ //if the var_name doesnot exist
    
	  index = mna_variable_indices.size();
    
	  mna_variable_indices[var_name]= index ;
	  
    }
    
    return index;
}

//this function adds a has map that maps every inductor to its added current, so that it would be easy when we use the mutual inductor
void Circuit::add_inductor_index(std::string inductor_name, double value){
    
    //check if the node is already added
    if(inductors.find(inductor_name)==inductors.end()){
	  int index;
	  inductor_name+= ".I";
	  index =  get_variable_index(inductor_name);
	  inductors[inductor_name]= std::pair<int,double>(index,value) ;
    }else{
	  std::string err = std::string("Element with name ")+inductor_name+ std::string(" already exists");
	  throw std::runtime_error(inductor_name);
    }
}

//returns node_index if exists, -1 if ground , -2 if not exists
int Circuit::get_variable_index(std::string var_name){
    int result = -2 ;
    if(var_name=="0" || var_name=="gnd"){
	result = -1; //means it is a ground
    }else{
	//search the map first, may be we can find it
	std::map< std::string , int>::iterator it=mna_variable_indices.find(var_name);
	
	if(it!=mna_variable_indices.end()){
	    result = it->second;
	}else{ //Another chance, may be it is in the alias map (a map that has nodes with the same indices)
	      std::map< std::string , std::string>::iterator  similar_nodes_iter;
	      
	      for(similar_nodes_iter = similar_nodes.begin(); similar_nodes_iter!=similar_nodes.end(); similar_nodes_iter++){
		  if(similar_nodes_iter->first == var_name){
		      result = get_variable_index(similar_nodes_iter->second); //we have to get the index of the similar node
		  }else if(similar_nodes_iter->second == var_name){
		      result = get_variable_index(similar_nodes_iter->first); //we have to get the index of the similar node
		  }
	      }
	  
	}
	
    }
    

    return result;
}

void Circuit::operator << (Element* E){
    	components.push_back(E);
	if(!E->is_linear()){
		is_linear = false;
	}
}
     
void Circuit::attach_elements(){
  
     //let us first sort the elements by their parsing order.
     components.sort([](Element* e1, Element* e2) { return e1->get_order_index() < e2->get_order_index(); });
     
     std::list<Element*>::iterator iter;

     
      
     //check if there is any component that needs to add an extra variable
     for(iter=components.begin(); iter!=components.end(); iter++){
	 
	  //first add the nodes of this element
	  (*iter)->add_my_nodes(this); //Note that components is a list, sho if an element appended elemets to the circuit (subcircuit) 
					//that wont affect the iterator. This is only for lists but others like vectors the iterator could change if we append to it
	  
     }
     
     int number_of_mna_variables = (int)mna_variable_indices.size();

     //now create the mna matrices
     G.create(number_of_mna_variables, number_of_mna_variables);
     C.create(number_of_mna_variables, number_of_mna_variables);
     J.create(number_of_mna_variables, number_of_mna_variables);
     B.create(number_of_mna_variables, 1);
     fx.create(number_of_mna_variables, 1);
     
     for(iter=components.begin(); iter!=components.end(); iter++){
	  (*iter)->write_stamp(G,C,this);  
     }
}




void Circuit::update_probes(double time, const double* solution){
    std::list<Probe*>::iterator iter;
    for(iter= Probes.begin(); iter!= Probes.end(); iter++){
	(*iter)->get_data(time, solution);
    }
}



///This function calls the source.update_B() to put the new values of the sources in the B vector at time 
void Circuit::update_sources(double time){
    std::list<Source*>::iterator iter;
    
    for(iter= sources.begin(); iter!= sources.end(); iter++){
	(*iter)->update_B( B , time); //this is a function of all sources to update the B vector
    }
}


///This function calls the NonLinElement.update_fx() to put the new values of the non linear expression in the fx vector at time 
///Note: that the function NonLinElement.update_fx() is only defined for nonlinear elements as they are also inherited from "NonLinElement" class
void Circuit::update_fx(const double* solution){
    std::list<NonLinElement*>::iterator iter;
    
    for(iter= Non_Linear_Elements.begin(); iter!= Non_Linear_Elements.end(); iter++){
	(*iter)->update_fx( fx , solution); //this is a function of all nonlinear elements to update the fx vector
    }
}

///This function calls the NonLinElement.update_J() to put the new values of the non linear Jacobian matrix in the G matrix at time 
///Note: that the function NonLinElement.update_J() is only defined for nonlinear elements as they are also inherited from "NonLinElement" class
void Circuit::update_J(const double* solution){
    std::list<NonLinElement*>::iterator iter;
    
    for(iter= Non_Linear_Elements.begin(); iter!= Non_Linear_Elements.end(); iter++){
	(*iter)->update_J( J , solution); 
    }
}

void Circuit::start_analysis(){
    attach_elements();
    
    std::vector<Analysis*>::iterator iter;
  
    
    for(iter= required_analysis.begin(); iter!= required_analysis.end(); iter++){
	//check if it is a DC solution to save it for other analysis
	if(DC* D = dynamic_cast<DC*> (*iter)){
	    DC_Analysis = D;
	    have_dc_solution = true;
	}
	//start simulating the required analysis
	(*iter)->simulate(G, C, J, B, fx, this);
    }
}

void Circuit::plot_probes(){
    std::list<Probe*>::iterator iter;
  
    for(iter= Probes.begin(); iter!= Probes.end(); iter++){
	(*iter)->plot();
    }
    
}

const BMatrix::Dense<double>& Circuit::get_dc_solution(){
  if(have_dc_solution){
	return DC_Analysis->get_solution();
  }else{
        DC_Analysis = new DC;
	DC_Analysis->simulate(G, C, J, B, fx, this);
	have_dc_solution = true;
	
	return DC_Analysis->get_solution();
  }
}
    
Element* Circuit::search_elements(std::string name){
    
    Element* results = nullptr;
    
    std::list<Element*>::iterator search = find_if(components.begin(), components.end(), [name](Element* el){return el->get_name()==name;} );
    if(search==components.end()){
      results = nullptr;
    }else{
      results = *search;
    }
    
    return results;
}
