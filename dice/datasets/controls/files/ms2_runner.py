import time
import os

start = time.time()
files = os.listdir(".")
counter = 0
for g in files:
    if g.endswith(".mzdb"):
        counter += 1
        print counter, ": extracting scans", g
        os.system("python ms2.py %s" % (g[:-5]))
stop = time.time()
print "Scan extraction of %d files in %.2f seconds." % (counter, stop - start)
