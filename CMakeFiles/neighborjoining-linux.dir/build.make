# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_SOURCE_DIR = /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining

# Include any dependencies generated for this target.
include CMakeFiles/neighborjoining-linux.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/neighborjoining-linux.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/neighborjoining-linux.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/neighborjoining-linux.dir/flags.make

CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o: CMakeFiles/neighborjoining-linux.dir/flags.make
CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o: src/neighborjoining.c
CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o: CMakeFiles/neighborjoining-linux.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o -MF CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o.d -o CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o -c /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/neighborjoining.c

CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/neighborjoining.c > CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.i

CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/neighborjoining.c -o CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.s

CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o: CMakeFiles/neighborjoining-linux.dir/flags.make
CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o: src/fasta_parser.c
CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o: CMakeFiles/neighborjoining-linux.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o -MF CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o.d -o CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o -c /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/fasta_parser.c

CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/fasta_parser.c > CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.i

CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/fasta_parser.c -o CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.s

CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o: CMakeFiles/neighborjoining-linux.dir/flags.make
CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o: src/binary_tree.c
CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o: CMakeFiles/neighborjoining-linux.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o -MF CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o.d -o CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o -c /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/binary_tree.c

CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/binary_tree.c > CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.i

CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/binary_tree.c -o CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.s

CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o: CMakeFiles/neighborjoining-linux.dir/flags.make
CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o: src/clusters_matrix.c
CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o: CMakeFiles/neighborjoining-linux.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o -MF CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o.d -o CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o -c /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/clusters_matrix.c

CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/clusters_matrix.c > CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.i

CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/clusters_matrix.c -o CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.s

CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o: CMakeFiles/neighborjoining-linux.dir/flags.make
CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o: src/cmdline.c
CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o: CMakeFiles/neighborjoining-linux.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o -MF CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o.d -o CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o -c /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/cmdline.c

CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/cmdline.c > CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.i

CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/src/cmdline.c -o CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.s

CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o: CMakeFiles/neighborjoining-linux.dir/flags.make
CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o: include/levenshtein/levenshtein.c
CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o: CMakeFiles/neighborjoining-linux.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o -MF CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o.d -o CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o -c /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/include/levenshtein/levenshtein.c

CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/include/levenshtein/levenshtein.c > CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.i

CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/include/levenshtein/levenshtein.c -o CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.s

# Object files for target neighborjoining-linux
neighborjoining__linux_OBJECTS = \
"CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o" \
"CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o" \
"CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o" \
"CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o" \
"CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o" \
"CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o"

# External object files for target neighborjoining-linux
neighborjoining__linux_EXTERNAL_OBJECTS =

neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/src/neighborjoining.c.o
neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/src/fasta_parser.c.o
neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/src/binary_tree.c.o
neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/src/clusters_matrix.c.o
neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/src/cmdline.c.o
neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/include/levenshtein/levenshtein.c.o
neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/build.make
neighborjoining-linux: CMakeFiles/neighborjoining-linux.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking C executable neighborjoining-linux"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/neighborjoining-linux.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/neighborjoining-linux.dir/build: neighborjoining-linux
.PHONY : CMakeFiles/neighborjoining-linux.dir/build

CMakeFiles/neighborjoining-linux.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/neighborjoining-linux.dir/cmake_clean.cmake
.PHONY : CMakeFiles/neighborjoining-linux.dir/clean

CMakeFiles/neighborjoining-linux.dir/depend:
	cd /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining /home/pavlos/research/PROJECTS/neighborjoining_repeated/neighborjoining/CMakeFiles/neighborjoining-linux.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/neighborjoining-linux.dir/depend

