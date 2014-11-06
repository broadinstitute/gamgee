#!/bin/bash

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/travis/gcc_4.9.1-1_amd64.deb
sudo apt-get remove cpp libffi-dev
sudo dpkg --install gcc_4.9.1-1_amd64.deb
