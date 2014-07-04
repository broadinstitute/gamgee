#!/bin/bash

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
sudo pip install cpp-coveralls
coveralls -b . -r . -e lib -e test -e testdata -t $1
