# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/moteki/C++/block-DDA

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/moteki/C++/block-DDA

# Include any dependencies generated for this target.
include mvp_fft/CMakeFiles/MBT_fft_init.dir/depend.make

# Include the progress variables for this target.
include mvp_fft/CMakeFiles/MBT_fft_init.dir/progress.make

# Include the compile flags for this target's objects.
include mvp_fft/CMakeFiles/MBT_fft_init.dir/flags.make

mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o: mvp_fft/CMakeFiles/MBT_fft_init.dir/flags.make
mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o: mvp_fft/MBT_fft_init.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/moteki/C++/block-DDA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o"
	cd /Users/moteki/C++/block-DDA/mvp_fft && g++-6   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o -c /Users/moteki/C++/block-DDA/mvp_fft/MBT_fft_init.cpp

mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.i"
	cd /Users/moteki/C++/block-DDA/mvp_fft && g++-6  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/moteki/C++/block-DDA/mvp_fft/MBT_fft_init.cpp > CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.i

mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.s"
	cd /Users/moteki/C++/block-DDA/mvp_fft && g++-6  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/moteki/C++/block-DDA/mvp_fft/MBT_fft_init.cpp -o CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.s

mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.requires:

.PHONY : mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.requires

mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.provides: mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.requires
	$(MAKE) -f mvp_fft/CMakeFiles/MBT_fft_init.dir/build.make mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.provides.build
.PHONY : mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.provides

mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.provides.build: mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o


# Object files for target MBT_fft_init
MBT_fft_init_OBJECTS = \
"CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o"

# External object files for target MBT_fft_init
MBT_fft_init_EXTERNAL_OBJECTS =

mvp_fft/libMBT_fft_init.a: mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o
mvp_fft/libMBT_fft_init.a: mvp_fft/CMakeFiles/MBT_fft_init.dir/build.make
mvp_fft/libMBT_fft_init.a: mvp_fft/CMakeFiles/MBT_fft_init.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/moteki/C++/block-DDA/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libMBT_fft_init.a"
	cd /Users/moteki/C++/block-DDA/mvp_fft && $(CMAKE_COMMAND) -P CMakeFiles/MBT_fft_init.dir/cmake_clean_target.cmake
	cd /Users/moteki/C++/block-DDA/mvp_fft && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MBT_fft_init.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
mvp_fft/CMakeFiles/MBT_fft_init.dir/build: mvp_fft/libMBT_fft_init.a

.PHONY : mvp_fft/CMakeFiles/MBT_fft_init.dir/build

mvp_fft/CMakeFiles/MBT_fft_init.dir/requires: mvp_fft/CMakeFiles/MBT_fft_init.dir/MBT_fft_init.cpp.o.requires

.PHONY : mvp_fft/CMakeFiles/MBT_fft_init.dir/requires

mvp_fft/CMakeFiles/MBT_fft_init.dir/clean:
	cd /Users/moteki/C++/block-DDA/mvp_fft && $(CMAKE_COMMAND) -P CMakeFiles/MBT_fft_init.dir/cmake_clean.cmake
.PHONY : mvp_fft/CMakeFiles/MBT_fft_init.dir/clean

mvp_fft/CMakeFiles/MBT_fft_init.dir/depend:
	cd /Users/moteki/C++/block-DDA && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/moteki/C++/block-DDA /Users/moteki/C++/block-DDA/mvp_fft /Users/moteki/C++/block-DDA /Users/moteki/C++/block-DDA/mvp_fft /Users/moteki/C++/block-DDA/mvp_fft/CMakeFiles/MBT_fft_init.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : mvp_fft/CMakeFiles/MBT_fft_init.dir/depend
