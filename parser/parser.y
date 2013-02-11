
%{
/*--------------This part is normal C++ code  that would appear before the Bison code-----*/
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <circuit.h>

#define YYDEBUG 1

extern "C" FILE *yyin;

void yyerror(const char *str)
{
        fprintf(stderr,"error: %s\n",str);
}
  
extern "C"
{
        int yyparse(void);
        int yylex(void);  
        int yywrap()
        {
                return 1;
        }
	int yylinno;
}

//declare the main crcuit class to add the elements to it
Circuit MainCircuit;

%}

/*------This part is the declaration of what are the types that the flex would return----------*/
%union{
	double dval;
	char* str;
}

%token <dval> NUMBER
%token <str> COMMENT
%token <str> NEWLINE


%%
/*--------------------This is the Bison code part--------------------*/
statments:
	|statments statment
	;
statment:
	| comment
	| number newline
	| newline
	;
comment: 
	| COMMENT NEWLINE{printf("Comment:%s\n",$1);}
	;

number: 
	| NUMBER
	| NUMBER NUMBER{
 			if(N==0){ //Iam in the first line
				N = $1;	
				if(N%2!=0) N=N+1; //add one vertex in case of odd number
				adjacency = new double[N*N];
				memset(adjacency,0,(N*N)*sizeof(double));
			}else{
				adjacency[N*(line_number-1)+(int)($1)] = $2;
				adjacency[N*(int)($1)+(line_number-1)] = $2; //this line to make sure it is symmetric
			}

			}
	;
newline: 
	| NEWLINE{
			if(N>0) line_number++;
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