import os
import sys
import shutil
from datetime import datetime


def time_stamp():
    return '[MutMap:{}]'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

def clean_cmd(cmd):
    return ' '.join(cmd.split())

def check_executables(commands, program_name='MutMap'):
    missing = []
    for command in commands:
        if command not in missing and shutil.which(command) is None:
            missing.append(command)

    if missing:
        sys.stderr.write(('{} Error: Missing required command(s): {}.\n').format(time_stamp(),
                                                                                  ', '.join(missing)))
        sys.stderr.write(('{} Please install them and ensure they are available in your PATH.\n').format(time_stamp()))
        sys.stderr.write(('{} {} has stopped before running the pipeline.\n').format(time_stamp(),
                                                                                     program_name))
        sys.exit(1)

def call_log(out_dir, name, cmd):
    print(time_stamp(), 
          '!!ERROR!! {}\n'.format(cmd), 
          flush=True)

    print('please check below:\n')

    log_path = '{}/log/{}.log'.format(out_dir, name)
    if not os.path.isfile(log_path):
        print('log file is not found: {}'.format(log_path))
        return

    with open(log_path) as log:
        for line in log:
            print(line, end='')
