#!/bin/bash

# Überprüfen und Ausgeben der BASE_PATH Variable
if [ -z "$BASE_PATH" ]; then
    echo "BASE_PATH ist im Unterskript nicht gesetzt"
else
    echo "BASE_PATH im Unterskript: $BASE_PATH"
fi

less "$BASE_PATH/bin/run_all.sh"
