import time

total_start = time.time()

from mz.API import mzdb
import sys
import os

sample_name = sys.argv[1]

raw = mzdb.File(sample_name + ".mzdb")
nist_ids = {}

if not os.path.exists(sample_name):
    os.mkdir(sample_name)
f = open(sample_name + "/scans.tab", 'w')

for (offset, (t, p)) in enumerate(raw.scan_list()):
    scan_num = offset + 1
    scan_type = raw.filters()[offset][1][:4]
    assert scan_type in ["ITMS", "FTMS"]
    pz = raw.extra_info(scan_num)["Charge State"]
    nist_id = ""
    peptide = ""
    nist_score = ""
    if scan_num in nist_ids:
        (nist_id, peptide, nist_score) = nist_ids[scan_num]
    print >> f, "%.3f\t%.4f\t%.0f\t%s\t%d\t%s\t%s\t%s" % (t, p, pz, scan_type, scan_num, nist_id, peptide, nist_score)
f.close()
raw.close()

total_stop = time.time()

print "Rawfile extracted in %.2f seconds" % (total_stop - total_start)
