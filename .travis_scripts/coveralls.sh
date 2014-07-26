#!/bin/bash

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
if [ "$CXX" == "clang++" ];
then
  sudo pip install cpp-coveralls
  coveralls -b . -r . -e lib -e test -e testdata -t ${COVERALLS_TOKEN}
fi
