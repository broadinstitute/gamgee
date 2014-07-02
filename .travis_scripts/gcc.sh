#!/bin/bash

sudo apt-get install -qq g++-4.8;
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 50;

