set(SIMULATOR_SOURCES Probes.cpp SourceFunc.cpp Source.cpp Analysis.cpp Inductor.cpp Capacitor.cpp Circuit.cpp resistor.cpp subcircuit.cpp plot.cc)

add_library(simulator STATIC ${SIMULATOR_SOURCES})

set(TEMPLATE_DIR "${CMAKE_SOURCE_DIR}/templates")
set_target_properties(simulator PROPERTIES COMPILE_FLAGS "-I\"${TEMPLATE_DIR}\"")