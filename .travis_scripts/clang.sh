#!/bin/bash

sudo apt-get install --allow-unauthenticated -qq clang-3.4 libstdc++-4.8-dev
sudo update-alternatives --install /usr/local/bin/clang++ clang++ /usr/bin/clang++-3.4 50
