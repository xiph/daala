#!/usr/bin/env python

from __future__ import print_function

import requests
import argparse
import os
import subprocess
import sys

if "check_output" not in dir( subprocess ): # duck punch it in!
    def f(*popenargs, **kwargs):
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = f

if 'DAALA_ROOT' not in os.environ:
    print("Please specify the DAALA_ROOT environment variable to use this tool.")
    sys.exit(1)

key = None
with open('secret_key','r') as keyfile:
    key = keyfile.read().strip()

if key is None:
    print("Could not open secret_key")
    sys.exit(1)

daala_root = os.environ['DAALA_ROOT']
os.chdir(daala_root)

parser = argparse.ArgumentParser(description='Submit test to arewecompressedyet.com')
parser.add_argument('-branch',default=None)
parser.add_argument('-prefix',default=None)
parser.add_argument('-master',action='store_true',default=False)
args = parser.parse_args()

if args.branch is None:
    args.branch = subprocess.check_output('git symbolic-ref -q --short HEAD',shell=True).strip()

if args.prefix is None:
    args.prefix = args.branch

commit = subprocess.check_output('git rev-parse HEAD',shell=True).strip()
short = subprocess.check_output('git rev-parse --short HEAD',shell=True).strip()
date = subprocess.check_output(['git','show','-s','--format=%ci',commit]).strip()
date_short = date.split()[0]
user = args.prefix
is_master = args.master

run_id = user+'-'+date_short+'-'+short

print('Creating run '+run_id)
r = requests.post("https://arewecompressedyet.com/submit/job", {'run_id': run_id, 'commit': commit, 'master': is_master, 'key': key})
print(r)
