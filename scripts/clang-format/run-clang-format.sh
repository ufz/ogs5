#!/bin/bash

# run this script at the source root directory
find . \( -name "*.h" -or -name "*.cpp" \) -print0 | xargs -0 clang-format -i -style=file
