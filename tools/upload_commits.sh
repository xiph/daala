#!/bin/bash
#
# USAGE: upload_commits.sh <num_commits> [<upload.py ARGS>]
# NOTE: -y and --send_mail are automatically passed to upload.py

if [ $# == 0 ]; then
  echo "usage: $0 <num_commits> [<upload.py ARGS>]"
  echo "NOTE: -y and --send_mail are automatically passed to upload.py"
  exit 1
fi

dir=$(dirname $0)

COMMITS=$1
shift

for i in $(seq $COMMITS -1 1)
do
    msg=$(git log --pretty=format:%s -n 1 HEAD~$((i - 1)))
    ${dir}/upload.py -y --send_mail $@ --rev HEAD~$i..HEAD~$((i - 1)) -m "$msg"
done
