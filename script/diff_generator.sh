#!/usr/bin/env bash
# Source https://github.com/Project-OSRM/osrm-backend

set -o errexit
set -o pipefail
set -o nounset

# Runs the Clang Formatter in parallel on the code base.
# Return codes:
#  - 1 there are files to be formatted
#  - 0 everything looks fine

# Get CPU count
OS=$(uname)
NPROC=1
if [[ $OS = "Linux" ]] ; then
    NPROC=$(nproc)
elif [[ ${OS} = "Darwin" ]] ; then
    NPROC=$(sysctl -n hw.physicalcpu)
fi

# macos does not have clang-foramt pre-installed
if [[ $OS = "Darwin" ]] ; then
    brew install clang-format
fi

# Discover clang-format
if type clang-format-9 2> /dev/null ; then
    # For linux the command is clang-format-9
    CLANG_FORMAT=clang-format-9
    V=$(clang-format --version)
    echo "clang-format is ${V}"
elif type clang-format 2> /dev/null ; then
    # For mac the command is clang-format
    CLANG_FORMAT=clang-format
    V=$(clang-format --version)
    echo "clang-format is ${V}"
else
    echo "No appropriate clang-format found"
    exit 1
fi

# Format cpp/h files inplace.
# If all files are well-formated, then there
# should be no changes.
find ./src -type f -name '*.h' -o -name '*.cpp' \
| xargs -I{} -P ${NPROC} ${CLANG_FORMAT} -i -style=file {}

# Format R files inplace.
# If all files are well-formated, then there
# should be no changes.
Rscript -e "install.packages(\"styler\");\
            library(\"styler\"); \
            style_dir(transformers=tidyverse_style(indent_by=4), exclude_files=\"./R/RcppExports.R\");" 
