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


#include "Inductor.h"
#include <math.h>

/*----------_Inductor Class -----------------*/
void Inductor::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){

    
    if(n1_index > -1){
	G.put(n1_index , current_index , 1);
	G.put(current_index , n1_index , -1);
    }
    
    if(n2_index > -1) {
        G.put(n2_index , current_index , -1);
	G.put(current_index , n2_index , 1);
    }
    
    C.put(current_index , current_index , value);
    
}


void Inductor::add_my_nodes(Circuit* circuit){
    n1_index = circuit->add_mna_variable(n1);
    n2_index = circuit->add_mna_variable(n2);
    
    //add extra variable for current
    std::string current = name + ".I";
    current_index = circuit->add_mna_variable(current);
    
    //we have to add the current element in case we need it for the mutual_inductances
    circuit->add_inductor_index(name,value);

  
    
}

/*----------_Mutual Inductor Class -----------------*/
void MututalInductor::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    //n1 and n2 are the inductors that have the mutual coupling between them
    
    std::pair<int,double> Inductor1 = circ->get_inductor_element_index(n1); //the index of the current of the first inductor
    std::pair<int,double> Inductor2 = circ->get_inductor_element_index(n2); //the index of the current of the second inductor
    
    double M = value/sqrt(Inductor1.second*Inductor2.second);
    C.put(Inductor1.first , Inductor2.first , -M);
    C.put(Inductor2.first , Inductor1.first , -M);
    
}


