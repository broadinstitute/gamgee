---
layout: post
title: "Setting up your development environment"
description: ""
category: setup
tags: [setup]
---
{% include JB/setup %}

### External dependencies

Use your package manager to install the latest versions of the following libraries:

  * [Boost Libraries](http://www.boost.org/)
  * [Cmake](http://www.cmake.org/)

If you are on a Mac, get [homebrew](http://brew.sh/) or [macports](http://www.macports.org/). Don't install these yourself, unless you really want to. 

#### Quick note for *homebrew* installed boost libraries

Homebrew creates all boost libraries in /usr/local/libboost_*.dylib and for the multi-threaded libraries (which we use) it creates a -mt extension to all of them (e.g. libboost_thread-mt.dylib). This will give you build errors. One solution is to remove the -mt of the names of all the boost libraries. Another is to create symlinks without the -mt to the ones with -mt. In particular boost_log will have both -mt and non -mt versions, we use the -mt one, so you will have to remove the non-mt one and rename or symlink the -mt one with libboost_log.dylib. We are aware of this issue and will eventually make a pull-request to homebrew's boost recipe to fix this situation for everyone.

[This is hopefully not an issue anymore, can someone check so we can remove it from the blog?]


### Building Gamgee

After cloning the repository, you are ready to build (we no longer use git submodules).

Gamgee uses the cmake build system in a very standard way. You can follow the usual cmake workflow to build your code:

	mkdir build      # create a directory for your build files
	cd build     
	cmake ..         # generate the makefile and check that your system has all the dependencies
	make             # compile the static library
	make gamgee_test # if you want to build the test binary


### Running Gamgee

The library will be produced inside your build directory on the following path:

    gamgee/libgamgee.a

Gamgee does not have an "install" target yet. We are working on it for the first public release.


### Testing Gamgee

Gamgee comes with an extensive suite of unit tests. To run all tests first you need to build the test binary with 

    $ make gamgee_test

The test binary is produced in `test/gamgee_test` inside your build directory.

You can now run your test manually **from the source root directory**
    
    $ cd ..
    $ build/test/gamgee_test

All the tests will try to access test files in the testdata directory, that is why you must be in the source root as your working directory. You can set up your IDE to use that as your working directory and run tests directly from Clion/Eclipse/Xcode.
