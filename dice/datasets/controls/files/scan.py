import time

total_start = time.time()

from mz.API import mzdb
import NoDB
import sys
import os
import cStringIO

sample_name = sys.argv[1]

raw = mzdb.File(sample_name + ".mzdb")

if not os.path.exists(sample_name):
    os.mkdir(sample_name)
NoDB.create(sample_name + "/scans")

for (offset, (t, p)) in enumerate(raw.scan_list()):
    scan_num = offset + 1
    scan_type = raw.filters()[offset][1][:4]
    assert scan_type in ["ITMS", "FTMS"]
    pz = raw.extra_info(scan_num)["Charge State"]
    o = cStringIO.StringIO()
    if scan_type == "FTMS":
        print >> o, "%.4f,%.0f,0" % (p, pz)
        for (mz, y, _s, z) in raw.cscan(raw.scanForTime(t)):
            print >> o, "%.4f,%.0f,%d" % (mz, y, z)
    else:
        print >> o, "%.4f,%.0f" % (p, pz)
        for (mz, y) in raw.scan(t):
            print >> o, "%.4f,%.0f" % (mz, y)
    NoDB.store(o.getvalue())
    o.close()
NoDB.close()
raw.close()

total_stop = time.time()

print "Rawfile extracted in %.2f seconds" % (total_stop - total_start)
output = open(sample_name + "_completed_scans.txt", 'w')
print >> output, "Rawfile extracted in %.2f seconds" % (total_stop - total_start)
output.close()
