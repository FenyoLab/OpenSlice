import time
import sys
import os

import time
import os

start = time.time()
files = os.listdir(".")
samples = [fname[:-5] for fname in files if fname.endswith(".mzdb")]
completed = [fname[:-20] for fname in os.listdir(".") if fname.endswith("_completed_scans.txt")]
todo_samples = set(samples).difference(completed)


# This level of indirection exists only because of SQLite's
# weird behaviour -- there is some weird memory leak or something
# so I defensively launch a process per file. Sad, but necessary.
def scan_launch(sample_name):
    os.system("python scan.py %s" % sample_name)


from multiprocessing import Pool

if __name__ == "__main__":
    p = Pool(9)
    p.map(scan_launch, todo_samples)
