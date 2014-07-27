
%{
#include <string>
#include <list>
#include <stdio.h>
#include "parser.hpp"


%}

%x incl

%%
\*.* {return COMMENT;}

[Vv] {return V_VARIABLE;}
\( {return LBRACKET;}
\) {return RBRACKET;}
\"  {return DOUBLE_QUOTE;}

^[ \t]*[rR][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return RESISTOR;}

^[ \t]*[lL][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return INDUCTOR;}

^[ \t]*[cC][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CAPACITOR;}

^[ \t]*[kK][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return MUTUALINDUCTOR;}

^[ \t]*[vV][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return VOLTAGESOURCE;}

^[ \t]*[iI][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CURRENTSOURCE;}

^[ \t]*[eE][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return E_ELEMENT;}

^[ \t]*[xX][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return SUBCKT_INSTANCE;}

\"(\\.|[^"])*\" {yylval.str=strdup(yytext); return QUOTED_STRING;}

[A-Za-z_]+[0-9]* {yylval.str=strdup(yytext); return STRING;}

^[ \t]*\.[tT][rR][aA][nN] {return TRAN;}

-?[0-9]+("."[0-9])*|-?"."[0-9]+|-?[0-9]+E[-+]?[0-9]+|-?[0-9]+"."[0-9]*E[-+]?[0-9]+|-?"."[0-9]+E[-+]?[0-9]+ { yylval.dval = atof(yytext) ; return DVALUE;} 

^[ \t]*\.[Ss][Uu][Bb][Cc][Kk][Tt] {return SUBCKT;}
^[ \t]*\.[eE][nN][dD][sS] {return END_SUBCKT;}

^[ \t]*\.[Pp][Rr][Ii][Nn][Tt][ \t]+[Tt][Rr][Aa][Nn] {return PRINT_TRAN;} /* Print statments*/


\n   {return NEWLINE;}                   /* ignore end of line */;
[ \t]+                  /* ignore whitespace */;



  /* START: Include Part */

^[ \t]*\.[Ii][Nn][Cc][Ll][Uu][Dd][Ee][ \t]+ BEGIN(incl);

<incl>[ \t]*      /* eat the whitespace */
<incl>[^ \t\n]+   { /* got the include file name */
            yyin = fopen( yytext, "r" );
     
                 if ( ! yyin ){
                     printf( "Can not open include file %s",yytext );
                     abort();
                 }
     
     		yypush_buffer_state(yy_create_buffer( yyin, YY_BUF_SIZE ));
     
                BEGIN(INITIAL);
}

<<EOF>> {
     	     //http://flex.sourceforge.net/manual/Multiple-Input-Buffers.html#Multiple-Input-Buffers
     	     //get the previoust parser from the stack. usefull when we have include file
     	     yypop_buffer_state();
     
             if ( !YY_CURRENT_BUFFER ){
                     yyterminate(); //terminate parsing
             }
	}
  /* END: Include Part */
%%