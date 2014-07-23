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

#ifndef CONTROLED_SOURCES_H
#define CONTROLED_SOURCES_H

#include "element.h""

/*-------------------The VCVS---------------*/
class VCVS : public FourTerminal
{
  private:
  double gain; //the gain of the VCVS
  int current_index;
  
  public:
    VCVS(std::string _inn1, std::string _inn2,std::string _out1, std::string _out2, double _gain): FourTerminal(_inn1,_inn2,_out1,_out2),gain(_gain){
	  name = std::string("E")+".+"+out1+".-"+out2;
    }
    
    VCVS(std::string _name, std::string _inn1, std::string _inn2,std::string _out1, std::string _out2, double _gain): FourTerminal(_inn1,_inn2,_out1,_out2),gain(_gain){
	  name = _name;
    }
      
    VCVS* clone(){return new VCVS(*this);};
    
    void write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ);
    bool is_linear(){return true;};
    
    ///add the node name to the circuit, 
    void add_my_nodes(Circuit* circuit); 
};

#endif // CONTROLED_SOURCES_H
