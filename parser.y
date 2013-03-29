
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


//declare the main crcuit class to add the elements to it
Circuit MainCircuit;

%}

%error-verbose
 
/*------This part is the declaration of what are the types that the flex would return----------*/
%union{
	double dval;
	char* str;
	int ival;
	std::list<char*>*  arg_list;
}

%token <dval> DVALUE
%token <str> STRING


%token <str> RESISTOR
%token <str> INDUCTOR
%token <str> CAPACITOR
%token <str> MUTUALINDUCTOR

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
	//| subcircuit_statment
	;
	
source:
	|voltagesource_statment
	;
	
command:
	| tran_statment
	;
	
resistor_statment:
	| RESISTOR node node DVALUE NEWLINE{
	      MainCircuit<< new resistor($1, $2, $3, $4);
	}
	;
	
inductor_statment:
	| INDUCTOR node node DVALUE NEWLINE{
	      MainCircuit<< new Inductor($1, $2, $3, $4);
	}
	;
	

capacitor_statment:
	| CAPACITOR node node DVALUE NEWLINE{
	      MainCircuit<< new Capacitor($1, $2, $3, $4);
	}
	;

voltagesource_statment:
	| VOLTAGESOURCE node node DVALUE{ //DC voltage source
	      MainCircuit<< new VoltageSource ($1, $2, $3 , new DCSource($4) ) ;
	}
	;
	
voltagesource_statment:
	| CURRENTSOURCE node node DVALUE{ //DC voltage source
	       MainCircuit<< new CurrentSource ($1, $2, $3 , new DCSource($4) ) ;
	}	
	;
	
tran_statment:
	| TRAN DVALUE DVALUE{
	      MainCircuit<< new transient(0, $3 , $2);
	}
	;
	
//subcircuit_statment:
//	|SUBCKT node_list{
	
//	}
//	;

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
	  $1->push_front($2);
	  $$ = $1;
	  printf("asasasas");
    }
    |node{
          $$ = new std::list<char*>;
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
	MainCircuit.start_analysis();
	MainCircuit.plot_probes();
	
	return 0;
} 