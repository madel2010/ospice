
%{
#include <stdio.h>
#include "parser.hpp"
%}



%%
%.* {yylval.str = yytext; return COMMENT;}
([0-9]+)|([0-9]*\.?[0-9]+)([eE][-+]?[0-9]+)? {yylval.dval =atof(yytext); return NUMBER;}
\n   {return NEWLINE;}                   /* ignore end of line */;
[ \t]+                  /* ignore whitespace */;



%%