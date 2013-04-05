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
     std::vector<Element*>::iterator it_el;
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

int Circuit::add_mna_variable(std::string node_name){
    
    if(node_name=="0" || node_name=="gnd"){
	return -1;
    }
    
     int index;
    //check if the node is already added
    if(mna_variable_indices.find(node_name)==mna_variable_indices.end()){
    
	  index = mna_variable_indices.size();
    
	  mna_variable_indices[node_name]= index ;
    }else{
	index = mna_variable_indices[node_name];
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

int Circuit::get_variable_index(std::string variable_name){
    if(variable_name=="0" || variable_name=="gnd"){
	return -1; //means it is a ground
    }else{
	return mna_variable_indices[variable_name];
    }
}

void Circuit::operator << (Element* E){
    	components.push_back(E);
	if(!E->is_linear()){
		is_linear = false;
	}
}
     
void Circuit::attach_elements(){
     std::vector<Element*>::iterator iter;

     
      
     //check if there is any component that needs to add an extra variable
     for(iter=components.begin(); iter!=components.end(); iter++){
	 
	  //Since this is the main circuit, we do not have to append any string to the nodes. Appending only done for the subcircuits. Check subcircuit.h
	  //Therefore we only need to create an empty vector with the size of the element terminsls and pass it to add_my_nodes_function
          std::vector<std::string> append_to_node_names ((*iter)->get_terminals_name().size());
	  
	  //first add the nodes of this element
	  (*iter)->add_my_nodes(this);
	  
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
    std::vector<Probe*>::iterator iter;
    for(iter= Probes.begin(); iter!= Probes.end(); iter++){
	(*iter)->get_data(time, solution);
    }
}



///This function calls the source.update_B() to put the new values of the sources in the B vector at time 
void Circuit::update_sources(double time){
    std::vector<Source*>::iterator iter;
    
    for(iter= sources.begin(); iter!= sources.end(); iter++){
	(*iter)->update_B( B , time); //this is a function of all sources to update the B vector
    }
}


///This function calls the NonLinElement.update_fx() to put the new values of the non linear expression in the fx vector at time 
///Note: that the function NonLinElement.update_fx() is only defined for nonlinear elements as they are also inherited from "NonLinElement" class
void Circuit::update_fx(const double* solution){
    std::vector<NonLinElement*>::iterator iter;
    
    for(iter= Non_Linear_Elements.begin(); iter!= Non_Linear_Elements.end(); iter++){
	(*iter)->update_fx( fx , solution); //this is a function of all nonlinear elements to update the fx vector
    }
}

///This function calls the NonLinElement.update_J() to put the new values of the non linear Jacobian matrix in the G matrix at time 
///Note: that the function NonLinElement.update_J() is only defined for nonlinear elements as they are also inherited from "NonLinElement" class
void Circuit::update_J(const double* solution){
    std::vector<NonLinElement*>::iterator iter;
    
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
    std::vector<Probe*>::iterator iter;
  
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
    
