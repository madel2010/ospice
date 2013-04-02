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

class SubCircuit : public Circuit
{
private:
	std::vector<std::string> terminals;
	std::string name;

public:
	SubCircuit(){

	}
};

#endif
