#!/bin/bash

# use Python-2.7 

Dir=`dirname $0`

python $Dir/filterByReadConfig.py "$@"

if [ "$?" != "0" ]; then
    echo "Python script failed!"
    exit 1
fi

#dos2unix $3

