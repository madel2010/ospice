# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/spider/C++_projects/My_Simulator

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/spider/C++_projects/My_Simulator

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/spider/C++_projects/My_Simulator/CMakeFiles /home/spider/C++_projects/My_Simulator/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/spider/C++_projects/My_Simulator/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named Ospice

# Build rule for target.
Ospice: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 Ospice
.PHONY : Ospice

# fast build rule for target.
Ospice/fast:
	$(MAKE) -f CMakeFiles/Ospice.dir/build.make CMakeFiles/Ospice.dir/build
.PHONY : Ospice/fast

#=============================================================================
# Target rules for targets named OspiceParser

# Build rule for target.
OspiceParser: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 OspiceParser
.PHONY : OspiceParser

# fast build rule for target.
OspiceParser/fast:
	$(MAKE) -f CMakeFiles/OspiceParser.dir/build.make CMakeFiles/OspiceParser.dir/build
.PHONY : OspiceParser/fast

#=============================================================================
# Target rules for targets named doc

# Build rule for target.
doc: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 doc
.PHONY : doc

# fast build rule for target.
doc/fast:
	$(MAKE) -f CMakeFiles/doc.dir/build.make CMakeFiles/doc.dir/build
.PHONY : doc/fast

#=============================================================================
# Target rules for targets named suitesparse

# Build rule for target.
suitesparse: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 suitesparse
.PHONY : suitesparse

# fast build rule for target.
suitesparse/fast:
	$(MAKE) -f CMakeFiles/suitesparse.dir/build.make CMakeFiles/suitesparse.dir/build
.PHONY : suitesparse/fast

#=============================================================================
# Target rules for targets named simulator

# Build rule for target.
simulator: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 simulator
.PHONY : simulator

# fast build rule for target.
simulator/fast:
	$(MAKE) -f src/CMakeFiles/simulator.dir/build.make src/CMakeFiles/simulator.dir/build
.PHONY : simulator/fast

#=============================================================================
# Target rules for targets named SymbolicC++

# Build rule for target.
SymbolicC++: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 SymbolicC++
.PHONY : SymbolicC++

# fast build rule for target.
SymbolicC++/fast:
	$(MAKE) -f SymbolicC++3-3.35/CMakeFiles/SymbolicC++.dir/build.make SymbolicC++3-3.35/CMakeFiles/SymbolicC++.dir/build
.PHONY : SymbolicC++/fast

lex.o: lex.cc.o
.PHONY : lex.o

# target to build an object file
lex.cc.o:
	$(MAKE) -f CMakeFiles/OspiceParser.dir/build.make CMakeFiles/OspiceParser.dir/lex.cc.o
.PHONY : lex.cc.o

lex.i: lex.cc.i
.PHONY : lex.i

# target to preprocess a source file
lex.cc.i:
	$(MAKE) -f CMakeFiles/OspiceParser.dir/build.make CMakeFiles/OspiceParser.dir/lex.cc.i
.PHONY : lex.cc.i

lex.s: lex.cc.s
.PHONY : lex.s

# target to generate assembly for a file
lex.cc.s:
	$(MAKE) -f CMakeFiles/OspiceParser.dir/build.make CMakeFiles/OspiceParser.dir/lex.cc.s
.PHONY : lex.cc.s

main.o: main.cpp.o
.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/Ospice.dir/build.make CMakeFiles/Ospice.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i
.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/Ospice.dir/build.make CMakeFiles/Ospice.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s
.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/Ospice.dir/build.make CMakeFiles/Ospice.dir/main.cpp.s
.PHONY : main.cpp.s

parser.o: parser.cpp.o
.PHONY : parser.o

# target to build an object file
parser.cpp.o:
	$(MAKE) -f CMakeFiles/OspiceParser.dir/build.make CMakeFiles/OspiceParser.dir/parser.cpp.o
.PHONY : parser.cpp.o

parser.i: parser.cpp.i
.PHONY : parser.i

# target to preprocess a source file
parser.cpp.i:
	$(MAKE) -f CMakeFiles/OspiceParser.dir/build.make CMakeFiles/OspiceParser.dir/parser.cpp.i
.PHONY : parser.cpp.i

parser.s: parser.cpp.s
.PHONY : parser.s

# target to generate assembly for a file
parser.cpp.s:
	$(MAKE) -f CMakeFiles/OspiceParser.dir/build.make CMakeFiles/OspiceParser.dir/parser.cpp.s
.PHONY : parser.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... Ospice"
	@echo "... OspiceParser"
	@echo "... doc"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... suitesparse"
	@echo "... simulator"
	@echo "... SymbolicC++"
	@echo "... lex.o"
	@echo "... lex.i"
	@echo "... lex.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... parser.o"
	@echo "... parser.i"
	@echo "... parser.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
