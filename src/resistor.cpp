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


#include "resistor.h"
#include "Sparse.h"
#include "MatrixBase.h"

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

using boost::lexical_cast;
using boost::bad_lexical_cast;
    
/*-------------Linear Resisor--------------*/
void resistor::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    
    if(n1_index > -1){
	G.add_to_entry(n1_index , n1_index , 1/value);
    }
    
    if(n2_index > -1) {
        G.add_to_entry(n2_index , n2_index , 1/value);
    }
    
    if(n2_index > -1 && n1_index >-1 ) {
	G.add_to_entry(n1_index , n2_index , -1/value);
	G.add_to_entry(n2_index , n1_index , -1/value);
    }
}


void resistor::add_my_nodes(Circuit* circuit){
    n1_index = circuit->add_mna_variable(n1);
    n2_index = circuit->add_mna_variable(n2);
}



/*-------------Non Linear Resisor--------------*/
nonlin_resistor::nonlin_resistor(std::string _n1, std::string _n2, const char* _curr_expression):TwoTerminal(_n1,_n2){

	     name = std::string("nlR")+".+"+n1+".-"+n2;

             //Remove all spaces from the expression to make it more compact when parsing it	     
	     Expression =  boost::regex_replace(std::string(_curr_expression),boost::regex("\\s+"),"");
	     
	     //add paranthesis to the expression. It is required for the  shunting_yard algorithm to make sure the expression is done.

             Expression = std::string("(") + Expression;
	     Expression += ")";
	     
	}
	
nonlin_resistor::nonlin_resistor(std::string _name, std::string _n1, std::string _n2, const char* _curr_expression):TwoTerminal(_n1,_n2){
	
	     name =_name;
	     
	     Expression = boost::regex_replace(std::string(_curr_expression),boost::regex("\\s+"),"");
	     
	     //add paranthesis to the expression. It is required for the  shunting_yard algorithm to make sure the expression is done
	     Expression = std::string("(") + Expression;
	     Expression += ")";
	     
	};

void nonlin_resistor::write_stamp(BMatrix::Sparse<double> &G, BMatrix::Sparse<double> &C, Circuit* circ){
    
    circ->add_NonLinElement(this);
    
    //call the shunt_yard algorithm to parse the nonlinear expression
    shunting_yard(circ);
    
    //Find dF/dx. We assume that all the tree nodes has been evaluated in the shunting_yard algorithm, and F has been evaluated
    std::vector< std::pair<Symbolic , int> >::iterator iter;
    for(iter=expression_depends_on.begin(); iter!=expression_depends_on.end(); iter++){
	dFdx.push_back(df(F,iter->first)); //deffierentaite F with respect to each variable it depends on in expression_depends_on
    }
}


void nonlin_resistor::add_my_nodes(Circuit* circuit ){
    n1_index = circuit->add_mna_variable(n1);
    n2_index = circuit->add_mna_variable(n2);
}


//The Shunting Yard algorithm is used to parse Mathematical Expressions
//This function is a modefied version of the code found in "http://en.wikipedia.org/wiki/Shunting_yard_algorithm"
      
//#define is_operator(c)  (c.compare("+") || c.compare("-") || c.compare("/") || c.compare("*") || c.compare("!") || c.compare("%") || c.compare("="))
#define is_operator(c) (c == '+' || c == '-' || c == '/' || c == '*' || c == '=' || c == '^')
#define is_function(c)  (strcasecmp(c.c_str(),"exp")==0 || strcasecmp(c.c_str(),"ln")==0 || strcasecmp(c.c_str(),"sin")==0 || strcasecmp(c.c_str(),"cos")==0 )


bool nonlin_resistor::is_numeric(const std::string& c){
    bool result;
    
    //First try to check if it is a number. if not then it might be a reference to node voltage or branch current
    try
    {
        lexical_cast<double>(c);
	result = true;
    }
    catch(boost::bad_lexical_cast& e)
    {
	result = false;
    }
    
    return result;
}

bool nonlin_resistor::is_voltage_current(const std::string& c , std::string& dependent){
 
    //if it is a reference to node, then this function returns the string dependent which is the name of the node that the expression depends on. Otherwise return dependet=""
    bool result;
    
    dependent = "";

    //Check if it is a node voltage or current
    boost::regex expr( "V\\(.+\\)|I\\(.+\\)" , boost::regex::icase) ;
    result = boost::regex_match( c, expr ) ;

    //if it is a reference to node voltage or branch current, then determine the name of the node
    if(result){
	//user forgot to define the node name
	if(c.size()==3){
		throw std::runtime_error("No node or branch is defined");
	}

        int pos = 2; //start after "V(" or "I("
	while( pos < c.size()-1 ){ //ind before closing ")"
		dependent+= c[pos];
		pos++;
	}
    }

    return result;
}

std::string nonlin_resistor::get_expression_token(int& pos){
    std::string token = "";
    
    //used to check if the first chr is a digit. If it is a digit, then all other chrs should also be digit
    bool first_chr_is_digit = false;
    
    //if it is "-" then check the char before it. if it is operator or "(" then it is a negative sign
    if(Expression[pos]=='-' && (Expression[pos-1]=='(' || is_operator(Expression[pos-1]) ) ){ 
	    //get the number becuase it is a negative sign
	    token = "-";
	    pos++;
	    while(pos < Expression.size() && !is_operator(Expression[pos]) && !is_function(token) && Expression[pos]!=')' ){
	    	token+= Expression[pos];
	    	pos++;
	    }
    }else if(is_operator(Expression[pos])){ //it is an operator
	    token = Expression[pos];
	    pos++;
    }else if(Expression[pos]=='('){ //it is a paranthesis
	    token = Expression[pos];
	    pos++;
    }else if(Expression[pos]==')'){ //it is a paranthesis
	    token = Expression[pos];
	    pos++;
    }else if(tolower(Expression[pos])=='v' || tolower(Expression[pos])=='i'){ //If it is a refernece to voltage node or branch current
	while(Expression[pos]!=')' ){
	    token+= Expression[pos];
	    pos++;
	    if(pos >= Expression.size()){
		std::string Err = "Wrong expression: ";
		Err+= Expression;
		throw std::runtime_error(Err);
	    } 
	}
	//Now add the closing bracket ")" because we did not consider it in the previous loop
	token += ")";
	pos++;
	
    }else{
	while(pos < Expression.size() && !is_operator(Expression[pos]) && !is_function(token) && Expression[pos]!=')' ){
	    token+= Expression[pos];
	    pos++;
	}
    }
    
    return token;
}
 
 // operators
// precedence   operators       associativity
// 4            !               right to left
// 3            * / %           left to right
// 2            + -             left to right
// 1            =               right to left
int nonlin_resistor::op_preced(const char* c)
{
    switch( (*c) )    {
        case '!':
            return 5;
        case '*':  case '/': case '%':
            return 4;
        case '+': case '-':
            return 3;
	case '^':
            return 2;
        case '=':
            return 1;
    }
    return 0;
}
 
bool nonlin_resistor::op_left_assoc(const char* c)
{
    switch( (*c) )    {
        // left to right
      case '*': case '/': case '%': case '+': case '-' : case '^':
            return true;
        // right to left
        case '=': case '!':
            return false;
    }
    return false;
}
 


Symbolic nonlin_resistor::do_operator( vector<Symbolic>::iterator left ,  vector<Symbolic>::iterator right , const char* op ){
   Symbolic result;
   
   if(*op=='+'){
      result = *left + *right;
   }else if(*op=='-'){
     result = *left - *right;
   }else if(*op=='/'){
     result = *left / *right;
   }else if(*op=='*'){
     result = *left * *right;
   }else if(*op=='^'){
     result = *left ^ *right;
   }
   
  
   return result;
   
}

Symbolic nonlin_resistor::do_operator( vector<Symbolic>::iterator left  , const char* op ){
   Symbolic result;
   
   if( strcasecmp(op,"exp")==0 ){
      result = exp(*left);
   }else if( strcasecmp(op,"sin")==0 ){
     result = sin(*left);
   }else if(strcasecmp(op,"cos")==0){
     result = cos(*left);
   }else if( strcasecmp(op,"ln")==0 ){
     result = ln(*left);
   }

   return result;
   
}

void nonlin_resistor::shunting_yard(Circuit* circ)
{
    int strpos = 0 , strend = Expression.size();
    
    //the output and operator stack
    //Note: we initialize it with the size of the Expression, becuase we know that we can not exceed the size of the Expression. 
    //Initializing the size provides better performance
    //std::vector<std::string> output(Expression.size()) ;
    std::vector<std::string> stack(Expression.size()); 
    
    std::vector<Symbolic> Tree; //this vector would hold the parsed expression tree
    Symbolic Temp_F;
    
    std::string c;
    std::string sc;          // used for record stack element
 
    //int output_counter = 0 ;
    int stack_counter = 0;
    
    //while the string is not finished
    while(strpos < strend)   {
      
	 //get a token from the expression.
	 //Note:: this function changes the strpos to point for the next charachter after the exctracted tocken
	 c = get_expression_token(strpos); 
	  

	  // If the token is a number (identifier), then add it to the output queue.
         if(is_numeric(c))  {
		 Temp_F = atof(c.c_str());
		 Tree.push_back(Temp_F);
	 }else if(is_voltage_current(c,sc)){
		 Temp_F = sc;
		 Tree.push_back(Temp_F);
		 //Next, we have to save the node names and references that the expression depends on, in order to use it latter.
		int node_index = circ->get_variable_index(sc);
		if(node_index<-1){ //no node found
		  std::string Err = std::string("Error in parsing nonlinear expression, node does not exist: ") + sc;
		  throw std::runtime_error(Err);
		}else{
		  expression_depends_on.push_back( std::pair<Symbolic , int>(Temp_F,node_index) );
		}
         }else if(is_function(c))   {
                stack[stack_counter] = c;
                stack_counter++;
         }else if(is_operator( (*c.c_str()) ))  {
                while(stack_counter > 0)    {
                    sc = stack[stack_counter - 1];
                    // While there is an operator token, op2, at the top of the stack
                    // op1 is left-associative and its precedence is less than or equal to that of op2,
                    // or op1 has precedence less than that of op2,
                    // Let + and ^ be right associative.
                    // Correct transformation from 1^2+3 is 12^3+
                    // The differing operator priority decides pop / push
                    // If 2 operators have equal priority then associativity decides.
                    if(is_operator( (*sc.c_str()) ) &&
                        ((op_left_assoc(c.c_str()) && (op_preced(c.c_str()) <= op_preced(sc.c_str()))) ||
                           (op_preced(c.c_str()) < op_preced(sc.c_str()))))   { 
                        // Pop op2 off the stack, onto the output queue;
                        stack_counter--;
		    
			//add to tree
			Temp_F = do_operator( Tree.end()-2 , Tree.end()-1, sc.c_str() );
			Tree.pop_back();
			Tree.pop_back();
			Tree.push_back(Temp_F);

                    }else{
                        break;
                    }
                }
                // push op1 onto the stack.
                stack[stack_counter] = c;
                stack_counter++;
            
	    // If the token is a left parenthesis, then push it onto the stack.
	    }else if(strcasecmp(c.c_str(),"(")==0)   {
                stack[stack_counter] = c;
                stack_counter++;
            // If the token is a right parenthesis:
	    }else if(strcasecmp(c.c_str(),")")==0){
                bool pe = false;
                // Until the token at the top of the stack is a left parenthesis,
                // pop operators off the stack onto the output queue
                while(stack_counter > 0)     {
                    sc = stack[stack_counter - 1];
                    if(strcasecmp(sc.c_str(),"(")==0 )    {
                        pe = true;
                        break;
                    }else{
                        stack_counter--;
			
			//add to tree
			Temp_F = do_operator( Tree.end()-2 , Tree.end()-1, sc.c_str() );
			Tree.pop_back();
			Tree.pop_back();
			Tree.push_back(Temp_F);
                    }
                }
                // If the stack runs out without finding a left parenthesis, then there are mismatched parentheses.
                if(!pe){
		      std::string Err = std::string("Error: parentheses mismatched in: ") + Expression;
		      throw std::runtime_error(Err);
                }
                // Pop the left parenthesis from the stack, but not onto the output queue.
                stack_counter--;
                // If the token at the top of the stack is a function token, pop it onto the output queue.
                if(stack_counter > 0){
                    sc = stack[stack_counter - 1];
                    if(is_function(sc)){
                         stack_counter--;
			
			 Temp_F = do_operator( Tree.end()-1 , sc.c_str() );
			 Tree.pop_back();
			 Tree.push_back(Temp_F);
                    }
                }
            }else{
                 std::string Err = std::string("Error in parsing nonlinear expression, Unknown token: ") + c;
		 throw std::runtime_error(Err);
            }

    }
    // When there are no more tokens to read:
    // While there are still operator tokens in the stack:
    while(stack_counter > 0)  {
        sc = stack[stack_counter - 1];
        
	if(strcasecmp(sc.c_str(),"(")==0 ||  strcasecmp(sc.c_str(),")")==0){
	      std::string Err = std::string("Error: parentheses mismatched in: ") + Expression;
	      throw std::runtime_error(Err);
        }
        stack_counter--;
    }
   
    //Make sure that we have only one expression at the end
    if(Tree.size()>1){
	std::string Err = "Something went wrong when parsing: ";
	Err+= Expression;
	throw std::runtime_error(Err);
    }else{
	F = Temp_F[0];
    }
}

void nonlin_resistor::update_fx(BMatrix::Dense<double>& fx, const double* solution){
      //firxt Evalue the expression at the specific time
      
      
      std::vector< std::pair<Symbolic , int> >::iterator iter;
      double value;
      Symbolic temp_F = F;

      for(iter=expression_depends_on.begin(); iter!=expression_depends_on.end(); iter++){
	  //get the value of the voltages that the expression depends on
	  value = solution[iter->second];
	  
	  temp_F = temp_F.subst(iter->first, value);
      }



      try
      {
	  value = lexical_cast<double>(temp_F);
      }
      catch(boost::bad_lexical_cast& e)
      {
	  throw std::runtime_error("Expression can not be converted to a value, something wrong");
      }
    
      if(n1_index > -1) fx.put(n1_index,0,value);
      if(n2_index > -1) fx.put(n2_index,0,-value);
}

//This function gives the result of J by updating the entries of the jacobian matrix 
void nonlin_resistor::update_J(BMatrix::Sparse<double>& J, const double* solution){
  
      std::vector< std::pair<Symbolic , int> >::iterator iter_depends_on;
      std::vector<Symbolic>::iterator iter_dFdx;
      
      double value;
      Symbolic temp;
      int dfdx_counter=0;
      
      for( iter_dFdx=dFdx.begin();  iter_dFdx!=dFdx.end(); iter_dFdx++){
	  
	  
	  temp = (*iter_dFdx);

	  for(iter_depends_on=expression_depends_on.begin(); iter_depends_on!=expression_depends_on.end(); iter_depends_on++){
	      temp = temp.subst(iter_depends_on->first , solution[iter_depends_on->second]);
	  }

	  try
	  {
	      value = lexical_cast<double>(temp);
	  }
	  catch(boost::bad_lexical_cast& e)
	  {
	      throw std::runtime_error("DF/DX can not be converted to a value, something wrong");
	  }
	  
	  if(n1_index > -1){
	      J.put(n1_index , expression_depends_on[dfdx_counter].second , value);
	  }
    
	  if(n2_index > -1) {
	      J.put(n2_index , expression_depends_on[dfdx_counter].second , value);
	  }
      
	  dfdx_counter++;
      }
      
}
