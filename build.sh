#!/bin/bash
set -e
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)"

pushd $DIR/neusomatic
  mkdir build
    pushd build
      cmake ..
      make
    popd
popd
