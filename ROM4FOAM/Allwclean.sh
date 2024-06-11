#!/bin/bash

# Ensure the script runs from its directory
cd "$(dirname "$0")" || exit 1

# Clean all directories in src/
cd src/ || exit 1

for d in */ ; do
    echo "Cleaning in directory: $d"
    cd "$d" || exit 1
    if [ -f "Make/files" ]; then
        wclean
        if [ $? -ne 0 ]; then
            echo "$d did not clean"
            exit 1
        fi
    fi
    cd .. || exit 1
done

cd .. || exit 1

# Clean all directories in applications/*/*
for d in applications/*/*/ ; do
    echo "Cleaning in directory: $d"
    cd "$d" || exit 1
    if [ -f "Make/files" ]; then
        wclean
        if [ $? -ne 0 ]; then
            echo "$d did not clean"
            exit 1
        fi
    fi
    cd - > /dev/null || exit 1
done

echo "Cleaning finished successfully"
