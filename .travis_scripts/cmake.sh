#!/bin/bash

wget https://www.broadinstitute.org/gatk/eng/travis/cmake_3.1.0-rc1-1_amd64.deb
sudo apt-get remove cmake cmake-data
sudo dpkg --install cmake_3.1.0-rc1-1_amd64.deb
