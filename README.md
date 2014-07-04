[![Build Status](https://travis-ci.org/broadinstitute/gamgee.svg?branch=master)](https://travis-ci.org/broadinstitute/gamgee)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gamgee/badge.png)](https://coveralls.io/r/broadinstitute/gamgee)


Gamgee
======

A C++ interface for key next generation sequencing file formats. This library
provides an API for handilng SAM, BAM, CRAM, VCF, BCF, BED, intervals and
interval_list files.

The implementation is currently incomplete, but is under development. The APIs
are free to change until we reach version 1. 

* Project website: http://broadinstitute.github.io/gamgee/
* Developer docs: http://broadinstitute.github.io/gamgee/doxygen/



Setting up development environment
----------------------------------
Use your package manager to install the latest versions of the following libraries: 

  * [Boost](http://www.boost.org/)

If you are on a Mac, get [homebrew](http://brew.sh/) or [macports](http://www.macports.org/). Don't install these yourself, unless you really want to. Homebrew separates [boost-build](https://github.com/Homebrew/homebrew/blob/master/Library/Formula/boost-build.rb) from [boost](https://github.com/Homebrew/homebrew/blob/master/Library/Formula/boost.rb), so make sure you install both.

#### Quick note for homebrew installed boost libraries

Homebrew creates all boost libraries in /usr/local/libboost_*.dylib and for the multi-threaded libraries (which we use) it creates a -mt extension to all of them (e.g. libboost_thread-mt.dylib). This will give you build errors. One solution is to remove the -mt of the names of all the boost libraries. Another is to create symlinks without the -mt to the ones with -mt. In particular boost_log will have both -mt and non -mt versions, we use the -mt one, so you will have to remove the non-mt one and rename or symlink the -mt one with libboost_log.dylib. We are aware of this issue and will eventually make a pull-request to homebrew's boost recipe to fix this situation for everyone.


Cloning the repository and updating it
--------------------------------------
Gamgee currently has the following dependencies as git submodules: 
  * [htslib](http://www.github.com/samtools/htslib)

For that reason, cloning Gamgee includes an extra step: 

    $ hub clone broadinstitute/gamgee (simple git clone also works)
    $ git submodule update --init --recursive

Every time you do a pull (or fetch/rebase) and notice that the *.gitmodule* file has changed, you need to update your submodules accordingly like so:

    $ git submodule update --recursive

Using git submodules to handle dependency management is currently an experiment. We may choose to move away from this structure in the near future, so stay tuned.



Building Gamgee
----------------
In any unix system and MacOS X Mavericks (10.9+) or newer, all you need is to run the boost-build system in the Foghorn source root:

    $ b2 

And *voila*. If you want more customization, you can configure your default compilers and options by creating a user-config.jam or site-config.jam file with them. For detailed instructions check the [Boost-build documentation](http://www.boost.org/boost-build2/doc/html/bbv2/overview/configuration.html)

In MacOS X older than Mavericks (10.9) you have to either configure your system to
use clang by default (by editing your user-config.jam or site-config.jam) or
specify it in the command line. 

    $ b2 toolset=clang cxxflags=-stdlib=libc++ linkflags=libc++

#### Building Tips

* Boost-build works very well with parallel compilation. Simply add -j N where N is the
number of threads you want to use to compile.
* *b2 clean* doesn't work the same way as *ant* or *mvn* clean. If you want to rebuild the whole project from scratch, use the -a option as in : 

    ```
$ b2 -a
    ```

* The default build is a debug build with no optimizations and tracing/debug information so you can
  debug your code. If you want to produce a release binary with all optimizations use: 

    ```
$ b2 variant=release
    ```

* Setup a user-config.jam file in your home directory (or in the Foghorn source root) with your
  compiler configurations so you don't have to see warnings or use the shell script "b" every time.
  Here is my example user-config.jam for MacOS and Unix: 

    ```
# MacOS
$ cat ~/user-config.jam
using clang : : : <cxxflags>-stdlib=libc++ <linkflags>-stdlib=libc++ ;
    ```

    ```
# Unix
$ cat ~/user-config.jam
using gcc ;
    ```


Running Gamgee
---------------
The binary is produced in *bin/&lt;toolset&gt;/&lt;variant&gt;/foghorn*, where toolset will be the compiler used and variant is either debug or release. Here are two examples: 

    bin/clang-darwin-4.2.1/debug/foghorn
    bin/gcc-4.9.0/release/foghorn


Testing Gamgee
---------------
Gamgee comes with an extensive suite of unit tests. To run all tests simply do:

    $ b2 test

The test binary is produced in *test/bin/&lt;toolset&gt;/&lt;variant&gt;/foghorn*, where toolset will be the compiler used and variant is either debug or release. Here are two examples:

    test/bin/clang-darwin-4.2.1/debug/run
    test/bin/gcc-4.9.0/release/run

You can run your test manually if you want to see any outputs from debug information you have inserted in your tests, or if you want to change the input files. The tests are run automatically every time you compile them with b2 test.

For more information about debugging tests and writing tests using the boost unit test framework,
please refer to the [documentation](http://www.boost.org/doc/libs/1_55_0/libs/test/doc/html/index.html)
