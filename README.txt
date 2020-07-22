README
------------------

The folder titled "proj-gl-files" refers to all of the necessary files for the OpenGL implementation. A brief excerpt from  the included main.cpp:

* This is simple testbed for my GLSL implementation.
 *
 * Usage: ./driver -i <input-file> [ -s <solution-file> ] [ -o <stats-file> ]
 *     <input-file>      File with commands to run
 *     <solution-file>   File with solution to compare with
 *     <stats-file>      Dump statistics to this file rather than stdout
 *
 * Only the -i is mandatory.  You must specify a test to run.  For example:
 *
 * ./driver -i 00.txt
 *
 * This will save the result to output.png.  You may compare this result with a
 * reference solution:
 *
 * ./driver -i 00.txt -s 00.png 

Running "./grading-script.py ./" will give a comprehensive diagnosis of how many of the test image files were successfuly rendered using the implementation. At the moment of this writing, it passes with 100% success.


-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------


The folder titled "proj-rt-files" contains all of the files for the Ray Tracer implementation. To run the tests, run ./ray_tracer -i <test-file>, where test-file is one of the provided test files. A brief excerpt from  the included main.cpp:

Usage: ./ray_tracer -i <test-file> [ -s <solution-file> ] [ -o <stats-file> ] [ -x <debug-x-coord> -y <debug-y-coord> ]

  Examples:

  ./ray_tracer -i 00.txt

  Renders the scene described by 00.txt.  Dumps the result to output.png.

  ./ray_tracer -i 00.txt -s 00.png

  Renders the scene described by 00.txt.  Dumps the result to output.png.
  Compares this to the solution in 00.png and dumps the error to diff.png.
  Outputs a measure of the error. 
