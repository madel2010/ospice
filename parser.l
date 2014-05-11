
%{
#include <string>
#include <list>
#include <stdio.h>
#include "parser.hpp"


%}



%%
\*.* {return COMMENT;}

^[rR][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return RESISTOR;}

^[lL][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return INDUCTOR;}

^[cC][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CAPACITOR;}

^[kK][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return MUTUALINDUCTOR;}

^[vV][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return VOLTAGESOURCE;}

^[iI][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CURRENTSOURCE;}

[A-Za-z_]+[0-9]* {yylval.str=strdup(yytext); return STRING;}

^\.[tT][rR][aA][nN] {return TRAN;}

-?[0-9]+("."[0-9])*|-?"."[0-9]+|-?[0-9]+E[-+]?[0-9]+|-?[0-9]+"."[0-9]*E[-+]?[0-9]+|-?"."[0-9]+E[-+]?[0-9]+ { yylval.dval = atof(yytext) ; return DVALUE;} 


\n   {return NEWLINE;}                   /* ignore end of line */;
[ \t]+                  /* ignore whitespace */;



%%