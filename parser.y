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

#include "Source.h"
#include "SourceFunc.h"
#include "Analysis.h"
#include "Probes.h"
#include <list>
#include <string>


//declare the list of crcuit classes to add the elements to it
//this is also usefull when we add subcircuits. The first element in this list is always the main circtui
std::list<Circuit*> Circuit_lists(1,new Circuit);
Circuit* CurrentCircuit = Circuit_lists.front();


#define YYDEBUG 1

extern "C" FILE *yyin;

void yyerror(const char *str)
{
        fprintf(stderr,"error: %s\n",str);
}
  
int yylex(void);

extern "C"
{
        int yyparse(void);
        int yywrap()
        {
                return 1;
        }
	int yylinno;
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
	| subcircuit_statment
	;
	
source:
	|voltagesource_statment
	;
	
command:
	| tran_statment
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
	| VOLTAGESOURCE node node DVALUE{ //DC voltage source
	      (*CurrentCircuit)<< new VoltageSource ($1, $2, $3 , new DCSource($4) ) ;
	}
	;
	
voltagesource_statment:
	| CURRENTSOURCE node node DVALUE{ //DC voltage source
	       (*CurrentCircuit)<< new CurrentSource ($1, $2, $3 , new DCSource($4) ) ;
	}	
	;
	
tran_statment:
	| TRAN DVALUE DVALUE{
	      (*CurrentCircuit)<< new transient(0, $3 , $2);
	}
	;
	
subcircuit_statment:
	|SUBCKT STRING node_list{
	    std::vector<std::string> terminals{ std::make_move_iterator(std::begin($3)), std::make_move_iterator(std::end($3)) };
	    Circuit_lists.push_back(new SubCircuit S1($2 , terminals));
	    CurrentCircuit = Circuit_lists.back();
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
	  free($1);
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
	
	//Start simulating the circuit
	CurrentCircuit->start_analysis();
	CurrentCircuit->plot_probes();
	
	return 0;
} 