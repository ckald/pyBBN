import os
import re
import time
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Automatic test results collector')
    parser.add_argument('--test', required=True)
    parser.add_argument('--modified', default=None)
    parser.add_argument('--regexp', default=r'.+')
    parser.add_argument('--file', default='kawano_output.dat')
    args = parser.parse_args()

    regex = re.compile(args.regexp)

    if args.modified:
        # Time in seconds since epoch for time, in which logfile can be unmodified.
        modified_since = time.time() - float(args.modified)

    observables = []

    for dirpath, dirnames, files in os.walk(args.test):
        for folder in dirnames:
            if regex.search(folder):
                datafile = os.path.join(dirpath, folder, args.file)
                data = None
                try:
                    with open(datafile) as f:
                        if args.modified:
                            modified = os.stat(datafile).st_mtime
                            if modified < modified_since:
                                continue

                        for line in reversed(f.readlines()):
                            if 'Observables:' in line:
                                data = line.rstrip('\n')

                    if data:
                        observables.append("{}\t{}".format(data, os.path.join(dirpath, folder)))
                except Exception, e:
                    print(e)

    print("\n".join(observables))
