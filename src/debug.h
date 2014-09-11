#ifndef __DEBUG_H
#define __DEBUG_H

// Simple debug helpers

/*
1->Information about Analysis
2->Information about elements
3->Information about circuit
4->
5->
6->
7->
8->
9->Information about Matrices
*/

#ifdef DEBUG
#define _DD(l) if (debug_level >= (l))

extern char debug_level;
#else
#define _D if (0) 
#define _DD(l) if (0) 
#endif

#endif /* __DEBUG_H */
