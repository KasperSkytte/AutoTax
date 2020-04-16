#!/bin/bash
#script meant to be run through the autotax docker container image
set -o errexit -o pipefail -o nounset
if [ ${PWD##*/} == "test" ]
then
  echo "Error: tests must be run from the root of the AutoTax git repository, not the test/ subfolder"
  exit 1
fi
bats -j -t $((`nproc`-2)) /opt/autotax/test
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Time elapsed: $duration!"
