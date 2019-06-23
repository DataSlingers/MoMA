#!/bin/bash -ue

set -o errexit
set -o pipefail
set -o nounset

bash script/diff_generator.sh

MSG="The following files have been modified:"
SUCC_MSG="Everything looks fine."
dirty=$(git ls-files --modified)

if [[ $dirty ]]; then
    echo $MSG
    echo $dirty
    exit 1
else
    echo $SUCC_MSG
    exit 0
fi
