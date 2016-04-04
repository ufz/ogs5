#!/bin/bash

# run this script at the source root directory

# for CLANG <= 3.7
#find . \( -name "*.h" -or -name "*.cpp" \) -print0 | xargs -0 clang-format -i -style=file

# for CLANG > 3.7
find . \( -name "*.h" -or -name "*.cpp" \) -print0 | xargs -0 clang-format -i -sort-includes=false -style=file

