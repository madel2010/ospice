#ifndef DefinedMatrices_H
#define DefinedMatrices_H

#include "Sparse.h"


inline BMatrix::Sparse<double> EYE(int m){
   BMatrix::Sparse<double> result(m,m);
   for(int i=0; i<m; i++){
      result.put(i,i,1.0);
   }
   return result;
}


#endif