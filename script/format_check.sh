#!/bin/bash -ue

set -o errexit
set -o pipefail
set -o nounset

bash script/diff_generator.sh

MSG="The following files have been modified:"
dirty=$(git ls-files --modified)

if [[ $dirty ]]; then
    echo $MSG
    echo $dirty
    exit 1
else
    exit 0
fi

