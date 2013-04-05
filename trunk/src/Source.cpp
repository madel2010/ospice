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
#include "Source.h"

/*---Current Source-------*/
void CurrentSource::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    n1_index = circ->get_variable_index(n1);
    n2_index = circ->get_variable_index(n2);
    
    circ->add_source(this);
}

void CurrentSource::add_my_nodes(Circuit* circuit){
    n1_index = circuit->add_mna_variable(n1);
    n2_index = circuit->add_mna_variable(n2);
}

double CurrentSource::update_B(BMatrix::Dense<double> &B, double time){
    double value = Func->get_value(time);
    if(n1_index > -1) B.put(n1_index,0,value);
    if(n2_index > -1) B.put(n2_index,0,-value);
}


/*---Voltage Source-------*/
void VoltageSource::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    n1_index = circ->get_variable_index(n1);
    n2_index = circ->get_variable_index(n2);
    
    //get the index of the extra variable of the current
    std::string current = name + ".I";
    current_index = circ->get_variable_index(current);
    
    //write the G matrix stamp
    if(n1_index > -1){
	G.put(current_index , n1_index , 1);
	G.put(n1_index, current_index , -1);
    }
    
    if(n2_index > -1){
	G.put(current_index , n2_index , -1);
	G.put(n2_index, current_index , 1);
    }
    
    circ->add_source(this);
}

    
void VoltageSource::add_my_nodes(Circuit* circuit ){
    n1_index = circuit->add_mna_variable(n1);
    n2_index = circuit->add_mna_variable(n2);
    
    //add extra variable for current
    std::string current =  name + ".I";
    current_index = circuit->add_mna_variable(current);
    
    
}

double VoltageSource::update_B(BMatrix::Dense<double> &B, double time){
    double value = Func->get_value(time);
    
    B.put(current_index,0,-value);
}
