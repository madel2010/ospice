
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
\"  {return DOUBLE_QUOTE;}

^[ \t]*[rR][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return RESISTOR;}

^[ \t]*[lL][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return INDUCTOR;}

^[ \t]*[cC][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CAPACITOR;}

^[ \t]*[kK][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return MUTUALINDUCTOR;}

^[ \t]*[vV][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return VOLTAGESOURCE;}

^[ \t]*[iI][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return CURRENTSOURCE;}

^[ \t]*[eE][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return E_ELEMENT;}

^[ \t]*[xX][A-Za-z0-9_]+ {yylval.str=strdup(yytext); return SUBCKT_INSTANCE;}

(?:.*[a-zA-Z]+[0-9]*)(?:.*\/).* {printf("%s",yytext);yylval.str=strdup(yytext); return PATH;}

[A-Za-z_]+[0-9]* {yylval.str=strdup(yytext); return STRING;}

^[ \t]*\.[tT][rR][aA][nN] {return TRAN;}

-?[0-9]+("."[0-9])*|-?"."[0-9]+|-?[0-9]+E[-+]?[0-9]+|-?[0-9]+"."[0-9]*E[-+]?[0-9]+|-?"."[0-9]+E[-+]?[0-9]+ { yylval.dval = atof(yytext) ; return DVALUE;} 

^[ \t]*\.[Ss][Uu][Bb][Cc][Kk][Tt] {return SUBCKT;}
^[ \t]*\.[eE][nN][dD][sS] {return END_SUBCKT;}

^[ \t]*\.[Pp][Rr][Ii][Nn][Tt][ \t]+[Tt][Rr][Aa][Nn] {return PRINT_TRAN;} /* Print statments*/

^[ \t]*\.[Ii][Nn][Cc][Ll][Uu][Dd][Ee] {return INCLUDE;} /* include statments*/





\n   {return NEWLINE;}                   /* ignore end of line */;
[ \t]+                  /* ignore whitespace */;


%%