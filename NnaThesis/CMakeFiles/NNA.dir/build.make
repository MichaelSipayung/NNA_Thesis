# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/michael/Documents/nna/NNA_Thesis/NNA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/michael/Documents/nna/NNA_Thesis/NNA/NnaThesis

# Include any dependencies generated for this target.
include CMakeFiles/NNA.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/NNA.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/NNA.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/NNA.dir/flags.make

CMakeFiles/NNA.dir/NNA.cpp.o: CMakeFiles/NNA.dir/flags.make
CMakeFiles/NNA.dir/NNA.cpp.o: ../NNA.cpp
CMakeFiles/NNA.dir/NNA.cpp.o: CMakeFiles/NNA.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/michael/Documents/nna/NNA_Thesis/NNA/NnaThesis/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/NNA.dir/NNA.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/NNA.dir/NNA.cpp.o -MF CMakeFiles/NNA.dir/NNA.cpp.o.d -o CMakeFiles/NNA.dir/NNA.cpp.o -c /home/michael/Documents/nna/NNA_Thesis/NNA/NNA.cpp

CMakeFiles/NNA.dir/NNA.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NNA.dir/NNA.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/michael/Documents/nna/NNA_Thesis/NNA/NNA.cpp > CMakeFiles/NNA.dir/NNA.cpp.i

CMakeFiles/NNA.dir/NNA.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NNA.dir/NNA.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/michael/Documents/nna/NNA_Thesis/NNA/NNA.cpp -o CMakeFiles/NNA.dir/NNA.cpp.s

# Object files for target NNA
NNA_OBJECTS = \
"CMakeFiles/NNA.dir/NNA.cpp.o"

# External object files for target NNA
NNA_EXTERNAL_OBJECTS =

NNA: CMakeFiles/NNA.dir/NNA.cpp.o
NNA: CMakeFiles/NNA.dir/build.make
NNA: CMakeFiles/NNA.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/michael/Documents/nna/NNA_Thesis/NNA/NnaThesis/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable NNA"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NNA.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/NNA.dir/build: NNA
.PHONY : CMakeFiles/NNA.dir/build

CMakeFiles/NNA.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NNA.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NNA.dir/clean

CMakeFiles/NNA.dir/depend:
	cd /home/michael/Documents/nna/NNA_Thesis/NNA/NnaThesis && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/michael/Documents/nna/NNA_Thesis/NNA /home/michael/Documents/nna/NNA_Thesis/NNA /home/michael/Documents/nna/NNA_Thesis/NNA/NnaThesis /home/michael/Documents/nna/NNA_Thesis/NNA/NnaThesis /home/michael/Documents/nna/NNA_Thesis/NNA/NnaThesis/CMakeFiles/NNA.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NNA.dir/depend

