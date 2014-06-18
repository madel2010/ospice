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
#include <map>


//declare the list of crcuit classes to add the elements to it
//this is also usefull when we add subcircuits. The first element in this list is always the main circtui
std::list<Circuit*> Circuit_lists(1,new Circuit);
Circuit* CurrentCircuit = Circuit_lists.front();

//<Name of instance , std::pair< name of subcircuit , terminals>
std::map<std::string , std::pair<std::string , std::list<std::string> > > subckt_instances;

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
}

%}

%error-verbose
 
/*------This part is the declaration of what are the types that the flex would return----------*/
%union{
	double dval;
	char* str;
	int ival;
	std::list<std::string>*  arg_list;
}
%destructor { delete $$; } <arg_list> 


%token <dval> DVALUE
%token <str> STRING


%token <str> RESISTOR
%token <str> INDUCTOR
%token <str> CAPACITOR
%token <str> MUTUALINDUCTOR

%token <str> SUBCKT_INSTANCE
%token <str> SUBCKT
%token <str> END_SUBCKT

%token <str> VOLTAGESOURCE
%token <str> CURRENTSOURCE

%token <str> TRAN

%token <str> COMMENT
%token <str> NEWLINE

%type <str> node
%type <arg_list> node_list
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
	;
	
source:
	|voltagesource_statment
	;
	
command:
	| tran_statment
	| subcircuit_statment
	| end_subckt_statment
	;
	
resistor_statment:
	| RESISTOR node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new resistor($1, $2, $3, $4);
	}
	;
	
inductor_statment:
	| INDUCTOR node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new Inductor($1, $2, $3, $4);
	}
	;
	

capacitor_statment:
	| CAPACITOR node node DVALUE NEWLINE{
	      (*CurrentCircuit)<< new Capacitor($1, $2, $3, $4);
	}
	;

voltagesource_statment:
	| VOLTAGESOURCE node node DVALUE NEWLINE{ //DC voltage source
	      (*CurrentCircuit)<< new VoltageSource ($1, $2, $3 , new DCSource($4) ) ;
	}
	;
	
voltagesource_statment:
	| CURRENTSOURCE node node DVALUE NEWLINE{ //DC voltage source
	       (*CurrentCircuit)<< new CurrentSource ($1, $2, $3 , new DCSource($4) ) ;
	}	
	;
	
tran_statment:
	| TRAN DVALUE DVALUE{
	      (*CurrentCircuit)<< new transient(0, $3 , $2);
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

%%


/*--------------------This is the main part of the code--------------------*/



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

	      
	      
	//Start simulating the circuit
	CurrentCircuit->start_analysis();
	CurrentCircuit->plot_probes();
	
	return 0;
} 