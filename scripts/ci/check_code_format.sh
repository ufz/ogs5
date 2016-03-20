#!/bin/bash
#
# Modified from https://github.com/google/closure-library/tree/master/scripts/ci
#
# Script to determine if files in Pull Request are properly formatted.
# Exits with non 0 exit code if formatting is needed.

CLANG_FORMAT_DIFF=clang-format-diff-3.7

FILES_TO_CHECK=$(git diff --name-only master | grep -E "\.h|\.cpp$")

if [ -z "${FILES_TO_CHECK}" ]; then
  echo "No files to check for formatting."
  exit 0
fi

FORMAT_DIFF=$(git diff -U0 master -- ${FILES_TO_CHECK} |
              ${CLANG_FORMAT_DIFF} -p1 -style=file)

if [ -z "${FORMAT_DIFF}" ]; then
  # No formatting necessary.
  echo "All files in PR properly formatted."
  exit 0
else
  # Found diffs.
  echo "ERROR: Found formatting errors!"
  echo "${FORMAT_DIFF}"
  exit 1
fi
