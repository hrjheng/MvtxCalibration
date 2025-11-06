#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
source /opt/sphenix/core/bin/sphenix_setup.sh -n new

export MYINSTALL=/sphenix/u/hjheng/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

# print the environment - needed for debugging
# printenv

ddump -s -g -n 0 -p $1 $2 | ./thrsana -f $3 -i 25 -n 5 -p $4 -t 1
# ddump -s -g -n $N -p 2001 $FILENAME | ./thrsana -f $FEEID -i 25 -n 5 -p $PREFIX$FEEID -t 1

echo all done
