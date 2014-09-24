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

*if you need detailed instructions, read our [blog post](http://broadinstitute.github.io/gamgee/setup/20140924-setting-up-your-development-environment/) on how to set up your development environment*.

You'll need the latest version of [Boost](http://www.boost.org/) installed in your system.

    $ hub clone broadinstitute/gamgee         # simple git clone also works here
    $ git submodule update --init --recursive # so you can get the submodules we depend on

Every time you do a pull (or fetch/rebase) and notice that the *.gitmodule* file has changed, you need to update your submodules accordingly like so:

    $ git submodule update --recursive

You can now build simply by running:

    $ b2

The binary is produced in *bin/&lt;toolset&gt;/&lt;variant&gt;/gamgee*, where toolset will be the compiler used and variant is either debug or release. Here are two examples:

    bin/clang-darwin-4.2.1/debug/gamgee
    bin/gcc-4.9.0/release/gamgee


Gamgee comes with an extensive suite of unit tests. To run all tests simply do:

    $ b2 test

The test binary is produced in *test/bin/&lt;toolset&gt;/&lt;variant&gt;/gamgee*, where toolset will be the compiler used and variant is either debug or release. Here are two examples:

    test/bin/clang-darwin-4.2.1/debug/run
    test/bin/gcc-4.9.0/release/run

You can run your test manually if you want to see any outputs from debug information you have inserted in your tests, or if you want to change the input files. The tests are run automatically every time you compile them with b2 test.
