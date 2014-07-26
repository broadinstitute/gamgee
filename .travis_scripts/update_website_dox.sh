#!/bin/bash

if [ "$CXX" == "clang++" ] && [ "$TRAVIS_BRANCH" == "master" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ]; 
then 
  echo -e "Downloading latest Doxygen...";
  cd ${HOME};
  wget https://github.com/doxygen/doxygen/archive/Release_1_8_7.tar.gz -O $HOME/doxygen-1.8.7.tgz;
  tar xzf doxygen-1.8.7.tgz;
  cd $HOME/doxygen-Release_1_8_7;
  ./configure > /dev/null;
  make > /dev/null;
  sudo make install > /dev/null;
  cd ${HOME}/build/broadinstitute/gamgee;
  doxygen

  echo -e "Publishing doxygen...\n";
  git config --global user.email "travis@travis-ci.org";
  git config --global user.name "travis-ci";
  git clone --branch=gh-pages https://${GH_TOKEN}@github.com/broadinstitute/gamgee gh-pages;
  cd gh-pages;
  rm -rf doxygen/;
  mv ../dox/html doxygen/;
  git add doxygen/;
  git commit -am "Latest doxygen documentation on successful travis build ${TRAVIS_BUILD_NUMBER} auto-pushed";
  git push origin gh-pages 

  echo -e "Published doxygen.\n"
  
fi

