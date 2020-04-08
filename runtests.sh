#!/bin/bash
if [ ${PWD##*/} == "test" ]
then
  echo "Error: tests must be run from the root of the AutoTax repository, not the test/ subfolder"
  exit 1
fi
sudo docker run --rm -it -v "$(pwd):/autotax" autotax:1.1 bats /autotax/test