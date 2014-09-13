%debug

%{
/*--------------This part is normal C++ code  that would appear before the Bison code-----*/


#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "resistor.h"
#include "Inductor.h"
#include "Capacitor.h"
#include "ControledSources.h"

#include "Circuit.h"
#include "subcircuit.h"

#include "Source.h"
#include "SourceFunc.h"
#include "Analysis.h"
#include "Probes.h"

#include <iterator> 
#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <boost/algorithm/string/trim.hpp>
#include <memory>




//declare the list of crcuit classes to add the elements to it
//this is also usefull when we add subcircuits. The first element in this list is always the main circtui
std::list<Circuit*> Circuit_lists(1,new Circuit);
Circuit* CurrentCircuit = Circuit_lists.front();

//<Name of instance , std::pair< name of subcircuit , terminals>
std::unordered_map<std::string , std::pair<std::string , std::list<std::string> > > subckt_instances;

//function to search for subcircuit name
struct FindSubckt: public std::binary_function< Circuit*, std::string, bool > {
		    bool operator () ( const Circuit* element, const std::string name ) const {
		    const SubCircuit* sub_circuit = dynamic_cast<const SubCircuit*>(element);
		    return sub_circuit->get_name() == name;
		}
};

#define YYDEBUG 1

extern "C" FILE *yyin;



  
int yylex(void);

extern "C"
{
        int yyparse(void);
        int yywrap()
        {
                return 1;
        }
	//int yylineno;
	//char* yytext;
}

void yyerror(const char *str){
        extern int yylineno;  
	extern char *yytext;
        printf("Line %d: %s at %s\n", yylineno, str, yytext);
        abort();
}

%}

%error-verbose
 
/*------This part is the declaration of what are the types that the flex would return----------*/
%union{
	double dval;
	/*char* str;*/
	char* str;
	int ival;
	std::list<std::string>*  arg_list;
}
%destructor { delete $$; } <arg_list> 
%destructor { delete $$; } <str> 

%token <str> EQUAL
%token <dval> DVALUE
%token <str> STRING
%token <str> QUOTED_STRING
%token <str> BEGIN_SPACE
%token <str> CUR

%token <str> RESISTOR
%token <str> INDUCTOR
%token <str> CAPACITOR
%token <str> MUTUALINDUCTOR
%token <str> E_ELEMENT
%token <str> F_ELEMENT
%token <str> G_ELEMENT

%token <str> SIN

%token <str> SUBCKT_INSTANCE
%token <str> SUBCKT
%token <str> END_SUBCKT

%token <str> VOLTAGESOURCE
%token <str> CURRENTSOURCE

%token <str> OP
%token <str> TRAN
%token <str> PRINT_TRAN
%token <str> PARAM

%token <str> COMMENT
%token <str> NEWLINE
%token <str> LBRACKET
%token <str> RBRACKET
%token <str> DOUBLE_QUOTE
%token <str> SINGLE_QUOTE
%token <str> V_VARIABLE
%token <str> I_VARIABLE

%type <str> node
%type <arg_list> node_list
%type <arg_list> v_node_list
%%
/*--------------------This is the Bison code part--------------------*/
statments:
	|statments statment
	;
statment:
	| comment
	| element newline
	| source newline
	| command newline
	| newline
	;
comment: 
	| COMMENT NEWLINE
	;

newline:
	| NEWLINE
	;
	
element:
	| resistor_statment
	| inductor_statment
	| capacitor_statment
	| subckt_instance_statment
	| vcvs_statment
	| cccs_statment
	;
	
source:
	|voltagesource_statment
	|currentsource_statment
	;
	
command:
	| op_statment
	| tran_statment
	| subcircuit_statment
	| end_subckt_statment
	| print_statment
	| param_statment
	;
	
resistor_statment:
	| RESISTOR node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new resistor(boost::trim_copy(std::string($1)), $2, $3, $4);
	}
	| G_ELEMENT node node CUR EQUAL QUOTED_STRING {
	      std::string expression = $6;
	      //remove the quotes at the begining and end 
	      expression.erase(expression.begin());
	      expression.erase(expression.end()-1);
	      (*CurrentCircuit)<<  new nonlin_resistor(boost::trim_copy(std::string($1)) , $2 , $3, expression.c_str());
	}
	;
		
	
inductor_statment:
	| INDUCTOR node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new Inductor(boost::trim_copy(std::string($1)), $2, $3, $4);
	}
	;
	

capacitor_statment:
	| CAPACITOR node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new Capacitor(boost::trim_copy(std::string($1)), $2, $3, $4);
	}
	;
	
vcvs_statment:
	| E_ELEMENT node node node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new VCVS(boost::trim_copy(std::string($1)), $4, $5, $2, $3,$6);
	}
	;

cccs_statment:
	| F_ELEMENT STRING node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new CCCS(boost::trim_copy(std::string($1)), $2, $3, $4, $5);
	}
	;
	
vccs_statment:
	| G_ELEMENT node node node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new VCCS(boost::trim_copy(std::string($1)), $4, $5, $2, $3,$6);
	}
	;	

voltagesource_statment:
	| VOLTAGESOURCE node node DVALUE NEWLINE{ //DC voltage source
	      (*CurrentCircuit)<< new VoltageSource (boost::trim_copy(std::string($1)), $2, $3 , new DCSource($4) ) ;
	}
	| VOLTAGESOURCE node node SIN LBRACKET DVALUE DVALUE DVALUE DVALUE DVALUE DVALUE RBRACKET{ //Sin voltage source
	      (*CurrentCircuit)<< new VoltageSource (boost::trim_copy(std::string($1)), $2, $3 , new SinSource($6,$7,$8,$9,$10,$11) ) ;
	}
	;
	
currentsource_statment:
	| CURRENTSOURCE node node DVALUE NEWLINE{ //DC voltage source
	       (*CurrentCircuit)<< new CurrentSource (boost::trim_copy(std::string($1)), $2, $3 , new DCSource($4) ) ;
	}	
	;

op_statment:
	| OP{
	      (*CurrentCircuit)<< new DC;
	}
	;
	
tran_statment:
	| TRAN DVALUE DVALUE{
	      (*CurrentCircuit)<< new transient(0, $3 , $2);
	}
	;
	
param_statment:
	| PARAM STRING EQUAL DVALUE {
	  std::ostringstream temp;
	  temp << $4;
	  CurrentCircuit->add_parameter($2 , temp.str());
	}
	;
	
subckt_instance_statment:
	| SUBCKT_INSTANCE node_list STRING{
	      subckt_instances[$1] = std::pair<std::string , std::list<std::string> >( $3 , (*$2) ) ;
	}
	;
	
subcircuit_statment:
	|SUBCKT STRING node_list NEWLINE{
	    //Initialize vecor from list. Note we are using C++11 syntax.
	    //{}  calls what is called an std::initializer_list
	    std::vector<std::string> terminals {std::make_move_iterator($3->begin()), std::make_move_iterator($3->end())};
	    Circuit_lists.push_back(new SubCircuit($2 , terminals));
	    CurrentCircuit = Circuit_lists.back();
	}
	;
	
end_subckt_statment:
	|END_SUBCKT NEWLINE{
	   CurrentCircuit = Circuit_lists.front();
	}
	;

	
print_statment:
       |PRINT_TRAN v_node_list{
	  for(std::string _node : (*$2)){
	    std::string node_name = _node.substr (2,_node.length()-3);
	    boost::trim(node_name);
	    if(_node[0]=='V'){ //Voltage Probale
	    	(*CurrentCircuit)<< new VoltageProbe(_node , node_name , "0");
	    }else if(_node[0]=='I' ){ //current Probe
		(*CurrentCircuit)<< new CurrentProbe(_node , node_name);
            }
	  }
       }
       ;
node:
     | DVALUE{
	      std::ostringstream n;
	      n << $1;
	      std::string name = n.str();
	      $$ = strdup(name.c_str());
	      
     }
     | STRING{
	      $$ = $1;
     }
     ;
    
node_list:
    |node_list node{ 
	  $1->push_back($2);
	  $$ = $1;
	  
    }
    |node{
          $$ = new std::list<std::string>;
          $$->push_back($1);
    }
    ;
    
v_node_list:
    |v_node_list V_VARIABLE LBRACKET node RBRACKET{
          $1->push_back(std::string("V(")+$4+std::string(")"));
	  $$ = $1;
    }
    |V_VARIABLE LBRACKET node RBRACKET{
          $$ = new std::list<std::string>;
          $$->push_back(std::string("V(")+$3+std::string(")"));
    }
    |v_node_list I_VARIABLE LBRACKET node RBRACKET{
          $1->push_back(std::string("I(")+$4+std::string(")"));
	  $$ = $1;
    }
    |I_VARIABLE LBRACKET node RBRACKET{
          $$ = new std::list<std::string>;
          $$->push_back(std::string("I(")+$3+std::string(")"));
    }
    ;    

%%


/*--------------------This is the main part of the code--------------------*/


char debug_level = 9;
int main(int argc, char** argv)
{

int yydebug = 0;
	if(argc<2) throw std::runtime_error("No file to parse");
	
	FILE *file;
	file = fopen(argv[1],"r");
	if(!file) throw std::runtime_error("Can not open the file");
		
	yyin = file;
	printf("Parsing file %s\n",argv[1]);
        yyparse();
	
	///START: create the subcircuit instances from the map of saved instances
	for(auto instance : subckt_instances){ //this is a new C++11 syntax
	      //find the subckt. We start from second becuase first one is main circuit
	      std::list<Circuit*>::iterator result = find_if( ++Circuit_lists.begin() , Circuit_lists.end(), std::bind2nd( FindSubckt(), instance.second.first ) );
	      SubCircuit* casted_result;
	      if(result!=Circuit_lists.end()){
		  casted_result = dynamic_cast<SubCircuit*>(*result);
	      }else{ //can not find the subcircuit
		  throw std::runtime_error(std::string("Subcircuit ")+instance.second.first+std::string(" is not defined"));
	      }
	      //create the termination as a vector
	      std::vector<std::string> inst_termination {std::make_move_iterator(instance.second.second.begin()), std::make_move_iterator(instance.second.second.end())};
	      (*CurrentCircuit)<< casted_result->create_instance(instance.first , inst_termination);
	}     
	///END: create the subcircuit instances from the map of saved instances

	      
	//Configuration of the circuit
	CurrentCircuit->config.files_base_name = argv[1];
	CurrentCircuit->config.output_mode = 1; //0->screen , 1->file

	//Start simulating the circuit
	CurrentCircuit->start_analysis();
	CurrentCircuit->plot_probes();
	
	return 0;
} 