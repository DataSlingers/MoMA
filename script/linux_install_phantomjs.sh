#!/usr/bin/env bash

# Reference: 
# https://rstudio.github.io/shinytest/articles/ci.html

OS=$(uname)

if [[ $OS = Linux ]] ; then
   export PHANTOMJS_VERSION=2.1.1
   phantomjs version
   export PATH=$PWD/travis_phantomjs/phantomjs$PHANTOMJS_VERSIONlinuxx86_64/bin:$PATH
   hash r
   phantomjs version
   if [[ $(phantomjs version) != $PHANTOMJS_VERSION ]]; then rm rf $PWD/travis_phantomjs; mkdir p $PWD/travis_phantomjs; fi
   if [[ $(phantomjs version) != $PHANTOMJS_VERSION ]]; then wget https://github.com/Medium/phantomjs/releases/download/v$PHANTOMJS_VERSION/phantomjs$PHANTOMJS_VERSIONlinuxx86_64.tar.bz2 O $PWD/travis_phantomjs/phantomjs$PHANTOMJS_VERSIONlinuxx86_64.tar.bz2; fi
   if [[ $(phantomjs version) != $PHANTOMJS_VERSION ]]; then tar xvf $PWD/travis_phantomjs/phantomjs$PHANTOMJS_VERSIONlinuxx86_64.tar.bz2 C $PWD/travis_phantomjs; fi
   if [[ $(phantomjs version) != $PHANTOMJS_VERSION ]]; then hash r; fi
   phantomjs version
fi
