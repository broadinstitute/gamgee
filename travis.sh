#!/bin/bash 

# add user-config.jam so boost-build uses the right GCC
echo "using gcc : 4.8 ;" > ~/user-config.jam

# build htslib (and all external dependencies)
cd lib/htslib && make lib-static && cd ../../

# build boost-build
cd lib/b2 && ./bootstrap.sh && ./b2 install --prefix=../../boost-build && cd ../../

