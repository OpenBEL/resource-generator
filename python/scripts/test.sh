#!/usr/bin/env bash
#
# Executes nosetests.
#
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"/../../
. "$DIR"/env.sh
cd "$RG_PYTHON" || exit 1

NOSE_OPTS+="--with-coverage "
NOSE_OPTS+="--cover-package=reggie "
NOSE_OPTS+="--cover-erase "
NOSE_OPTS+="--nocapture "
NOSE_OPTS+="--with-doctest "
NOSE_OPTS+="--doctest-options=+ELLIPSIS "

PYTHONPATH=$(pwd)/tests nosetests $NOSE_OPTS -v