#!/bin/bash

# Ensure the script runs from its directory
cd "$(dirname "$0")" || exit 1

# Compile all directories in src/
cd src/ || exit 1

for d in */ ; do
    echo "Compiling in directory: $d"
    cd "$d" || exit 1
    if [ -f "Make/files" ]; then
        wmake "$1" "$2"
        if [ $? -ne 0 ]; then
            echo "$d did not compile"
            exit 1
        fi
    fi
    cd .. || exit 1
done

cd .. || exit 1

# Compile all directories in applications/*/*
for d in applications/*/*/ ; do
    echo "Compiling in directory: $d"
    cd "$d" || exit 1
    if [ -f "Make/files" ]; then
        wmake "$1" "$2"
        if [ $? -ne 0 ]; then
            echo "$d did not compile"
            exit 1
        fi
    fi
    cd - > /dev/null || exit 1
done

echo "Compilation finished successfully"
