[![Build Status](https://travis-ci.org/broadinstitute/gamgee.svg?branch=master)](https://travis-ci.org/broadinstitute/gamgee)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gamgee/badge.png)](https://coveralls.io/r/broadinstitute/gamgee)
[![Issue Stats](http://issuestats.com/github/broadinstitute/gamgee/badge/pr)](http://issuestats.com/github/broadinstitute/gamgee)
[![Issue Stats](http://issuestats.com/github/broadinstitute/gamgee/badge/issue)](http://issuestats.com/github/broadinstitute/gamgee)

Gamgee
======

A C++ interface for key next generation sequencing file formats. This library
provides an API for handilng SAM, BAM, CRAM, VCF, BCF, BED, intervals and
interval_list files.

The implementation is currently incomplete, but is under development. The APIs
are free to change until we reach version 1.

* Project website: http://broadinstitute.github.io/gamgee/
* Developer docs: http://broadinstitute.github.io/gamgee/doxygen/


Gamgee is licensed under the [MIT License](http://www.google.com.br/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0CB8QFjAA&url=http%3A%2F%2Fopensource.org%2Flicenses%2FMIT&ei=CPEiVNH-DoLxgwSG3oLwAg&usg=AFQjCNHLnploR8MB1GvOk6dHWdThFdiIbg&bvm=bv.75775273,d.eXY)



### Quick start

*If you need detailed instructions, read our [blog post](http://broadinstitute.github.io/gamgee/setup/20140924-setting-up-your-development-environment/) on how to set up your development environment*.

You'll need the latest version of [Boost](http://www.boost.org/) installed in your system.

    $ hub clone broadinstitute/gamgee         # simple git clone also works here
    
Building Gamgee
----------------
Gamgee is configured to use CMake for better integration with the [CLion](https://www.jetbrains.com/clion/) IDE. We recommend an out-of-source build. From your Gamgee root directory, create a build directory:

    $ mkdir build/

Navigate to that directory to complete the initial CMake setup:

    $ cd build/
    $ cmake ..

If building on the Broad servers, you will need extra arguments:

    $ cmake .. -DCMAKE_C_COMPILER=`which gcc` -DCMAKE_CXX_COMPILER=`which g++`

After the setup of the build configuration and makefiles, build Gamgee by executing make from the build/ directory:

    $ make

By default, CMake will execute the most recent explicitly specified build configuration (release/debug). To specify a build configuration use:

    $ make release

or

    $ make debug

for the release or debug configurations, respectively.

### Building Tips
 
* To clean your build, use 

    ```
    $ make clean
    ```
* If building on the Broad cluster you will need these dotkits
 
    ```
reuse -q .boost-1.55.0
reuse -q .fftw-3.3.4
reuse -q GCC-4.9
reuse -q CMake
    ```
NOTE: the CMake dotkit is only availble for systems running RHEL 6 or higher.  For other systems, you will have to compile your own version of CMake from the latest [source](http://www.cmake.org/download/)

Testing Gamgee
---------------
Gamgee comes with an extensive suite of unit tests. To run all tests simply do:

    $ make run_test

Setting up CLion
----------------
1. Download the latest version of the CLion IDE from [here](https://www.jetbrains.com/clion/)

2. Click through the walk-through and set up your preferences as desired.

3. At the startup menu, select "Import Project From Sources"

4. Navigate to the Gamgee repo directory and select it.

5. CLion will inform you that "Directory gamgee contains CMakeLists.txt".  Choose "Open Project".

6. CLion may not recognize all of the namespaces in the code until you execute a build from the IDE.  Select "Build" from the Run menu or click the build button in the upper right (green down arrow with ones and zeros).

7. If CLion still can't find all dependencies, go to File > Invalidate Caches / Restart...  This will restart CLion. Reindexing files and references may take a minute or two.

NOTE: Using CMake from CLion uses a different build directory than using CMake from the command line.
