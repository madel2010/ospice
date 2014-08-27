//The Shunting Yard algorithm is used to parse Mathematical Expressions
//This function is a modefied version of the code found in "http://en.wikipedia.org/wiki/Shunting_yard_algorithm"
      

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <Circuit.h>

//#define is_operator(c)  (c.compare("+") || c.compare("-") || c.compare("/") || c.compare("*") || c.compare("!") || c.compare("%") || c.compare("="))
#define is_operator(c) (c == '+' || c == '-' || c == '/' || c == '*' || c == '=' || c == '^')
#define is_function(c)  (strcasecmp(c.c_str(),"exp")==0 || strcasecmp(c.c_str(),"ln")==0 || strcasecmp(c.c_str(),"sin")==0 || strcasecmp(c.c_str(),"cos")==0 )


bool is_numeric(const std::string& c){
    bool result;
    
    //First try to check if it is a number. if not then it might be a reference to node voltage or branch current
    try
    {
        boost::lexical_cast<double>(c);
	result = true;
    }
    catch(boost::bad_lexical_cast& e)
    {
	result = false;
    }
    
    return result;
}

bool is_voltage_current(const std::string& c , std::string& dependent){
 
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

std::string get_expression_token(const std::string& Expression , int& pos){
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
int op_preced(const char* c)
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
 
bool op_left_assoc(const char* c)
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
 


Symbolic do_operator( vector<Symbolic>::iterator left ,  vector<Symbolic>::iterator right , const char* op ){
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

Symbolic do_operator( vector<Symbolic>::iterator left  , const char* op ){
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


Symbolic shunting_yard(const std::string& Expression , std::vector< std::pair<Symbolic , int> >& expression_depends_on , Circuit* circ)
{
    int strpos = 0 , strend = Expression.size();
    
    //the output and operator stack
    //Note: we initialize it with the size of the Expression, becuase we know that we can not exceed the size of the Expression. 
    //Initializing the size provides better performance
    //std::vector<std::string> output(Expression.size()) ;
    std::vector<std::string> stack(Expression.size()); 
    
    std::vector<Symbolic> Tree; //this vector would hold the parsed expression tree
    Symbolic F, Temp_F;
    
    std::string c;
    std::string sc;          // used for record stack element
 
    //int output_counter = 0 ;
    int stack_counter = 0;
    
    //while the string is not finished
    while(strpos < strend)   {
      
	 //get a token from the expression.
	 //Note:: this function changes the strpos to point for the next charachter after the exctracted tocken
	 c = get_expression_token(Expression , strpos); 
	  

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

    return F;
}
