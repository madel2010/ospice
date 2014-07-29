/*
 * <one line to give the library's name and an idea of what it does.>
 * Copyright 2014  <copyright holder> <email>
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License or (at your option) version 3 or any later version
 * accepted by the membership of KDE e.V. (or its successor approved
 * by the membership of KDE e.V.), which shall act as a proxy
 * defined in Section 14 of version 3 of the license.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#include "ControledSources.h"

/*-------------------The VCVS (E-element)---------------*/
void VCVS::add_my_nodes(Circuit* circuit){
    in1_index = circuit->add_mna_variable(in1);
    in2_index = circuit->add_mna_variable(in2);
    
    out1_index = circuit->add_mna_variable(out1);
    out2_index = circuit->add_mna_variable(out2);
    
    //add extra variable for current
    std::string current =  name + ".I";
    current_index = circuit->add_mna_variable(current);
}

void VCVS::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    if(out1_index > -1){
	G.add_to_entry(out1_index , current_index , 1);
	G.add_to_entry(current_index , out1_index , 1);
    }
    
    if(out2_index > -1) {
        G.add_to_entry(out2_index , current_index , -1);
	G.add_to_entry(current_index , out2_index , -1);
    }
    
    if(in1_index > -1){
	G.add_to_entry(current_index , in1_index , -gain);
    }
    
    if(in2_index > -1){
	G.add_to_entry(current_index , in2_index , gain);
    } 
}

/*-------------------The VCCS (G-element)---------------*/
void VCCS::add_my_nodes(Circuit* circuit){
    in1_index = circuit->add_mna_variable(in1);
    in2_index = circuit->add_mna_variable(in2);
    
    out1_index = circuit->add_mna_variable(out1);
    out2_index = circuit->add_mna_variable(out2);

}

void VCCS::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    if(out1_index > -1 && in1_index > -1 ){
	G.add_to_entry(out1_index , in1_index , gain);
    }
    
    if(out1_index > -1 && in2_index > -1 ){
	G.add_to_entry(out1_index , in2_index , -gain);
    }
    
    if(out2_index > -1 && in1_index > -1 ){
	G.add_to_entry(out2_index , in1_index , -gain);
    }
    
    if(out2_index > -1 && in2_index > -1 ){
	G.add_to_entry(out2_index , in1_index , gain);
    }
    
    
}
