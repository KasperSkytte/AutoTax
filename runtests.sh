#!/bin/bash
if [ ${PWD##*/} == "test" ]
then
  echo "Error: tests must be run from the root of the AutoTax git repository, not the test/ subfolder"
  exit 1
fi
sudo docker build -t kasperskytte/autotax:latest docker/
(sudo docker run --rm -it -v "$(pwd):/autotax" kasperskytte/autotax:latest bats -t -j $((`nproc`-2)) /autotax/test) |& tee test_result.log
duration=$(printf '%02dh:%02dm:%02ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Done in: $duration!"