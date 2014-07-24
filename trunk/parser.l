
%{
#include <string>
#include <list>
#include <stdio.h>
#include "parser.hpp"


%}



%%
\*.* {return COMMENT;}

[Vv] {return V_VARIABLE;}
\( {return LBRACKET;}
\) {return RBRACKET;}

[rR][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return RESISTOR;}

[lL][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return INDUCTOR;}

[cC][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CAPACITOR;}

[kK][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return MUTUALINDUCTOR;}

[vV][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return VOLTAGESOURCE;}

[iI][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CURRENTSOURCE;}

[eE][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return E_ELEMENT;}

[xX][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return SUBCKT_INSTANCE;}


[A-Za-z_]+[0-9]* {yylval.str=strdup(yytext); return STRING;}

\.[tT][rR][aA][nN] {return TRAN;}

-?[0-9]+("."[0-9])*|-?"."[0-9]+|-?[0-9]+E[-+]?[0-9]+|-?[0-9]+"."[0-9]*E[-+]?[0-9]+|-?"."[0-9]+E[-+]?[0-9]+ { yylval.dval = atof(yytext) ; return DVALUE;} 

\.[Ss][Uu][Bb][Cc][Kk][Tt] {return SUBCKT;}
\.[eE][nN][dD][sS] {return END_SUBCKT;}

\.[Pp][Rr][Ii][Nn][Tt][ \t]+[Tt][Rr][Aa][Nn] {return PRINT_TRAN;} /* Print statments*/





\n   {return NEWLINE;}                   /* ignore end of line */;
[ \t]+                  /* ignore whitespace */;


%%