set(SYMBOLICC++_SOURCES src/array.cpp src/cloning.cpp src/derive.cpp src/identity.cpp 
 src/matnorm.cpp src/matrix.cpp src/multinomial.cpp src/polynomial.cpp 
 src/quatern.cpp src/rational.cpp src/vecnorm.cpp src/vector.cpp 
 src/verylong.cpp src/constants.cpp src/equation.cpp src/functions.cpp 
 src/integrate.cpp src/number.cpp src/product.cpp src/solve.cpp src/sum.cpp 
 src/symbol.cpp src/symbolic.cpp src/symbolicc++.cpp src/symerror.cpp 
 src/symmatrix.cpp)

add_library(SymbolicC++ STATIC ${SYMBOLICC++_SOURCES})

set_target_properties(SymbolicC++ PROPERTIES COMPILE_FLAGS "-I\"${CMAKE_CURRENT_SOURCE_DIR}/include\" -fno-elide-constructors")

