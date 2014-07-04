#!/bin/bash

# install boost libraries
sudo apt-get install -qq boost1.55 

# install boost build from sources in ${HOME}/boost-build
git clone https://github.com/boostorg/build.git ${HOME}/boost-build
cd ${HOME}/boost-build
./bootstrap.sh --with-toolset=$1
sudo ./b2 install 
cd -
