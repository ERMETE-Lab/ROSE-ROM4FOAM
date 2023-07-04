#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

cd src/

for d in ./*/ ; do   
    echo $d    
    cd $d
    if [ -f "Make/files" ]; then
    wclean $1 $2
    fi
    cd ../
done

cd ..

for d in applications/*/* ; do   
    echo $d    
    cd $d
    if [ -f "Make/files" ]; then
    wclean $1 $2
    if [ $? -ne 0 ]
    then
     echo $d "did not compile"
     exit 1
    fi
	fi
	cd ../../../
	pwd
done
