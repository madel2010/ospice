#ifndef _SHUNTING_YARD_H
#define _SHUNTING_YARD_H

#include "symbolicc++.h"

//The Shunting Yard algorithm is used to parse Mathematical Expressions
//This function is a modefied version of the code found in "http://en.wikipedia.org/wiki/Shunting_yard_algorithm"
Symbolic shunting_yard(const std::string& Expression , std::vector< std::pair<Symbolic , int> >& expression_depends_on , Circuit* circ);
      
//the following functions are used for shunting Yard algorithm
bool is_numeric(const std::string& c); //if it is a numeric token
bool is_voltage_current(const std::string& c , std::string& dependent); //check if it is a voltage or current token
std::string get_expression_token(const std::string& Expression, int& pos);

int op_preced(const char* c);

bool op_left_assoc(const char* c);

Symbolic do_operator( vector<Symbolic>::iterator left ,  vector<Symbolic>::iterator right , const char* op );

Symbolic do_operator( vector<Symbolic>::iterator left , const char* op );

#endif