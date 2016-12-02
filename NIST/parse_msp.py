import NoDB
import cStringIO
from sys import argv

msp_file = argv[1]
f = open(msp_file)
NoDB.create(msp_file[:-4])
entry_num = 0
nist_id = 0
while True:
    l = f.readline()
    if not l:
        break
    if l.startswith("ID: "):
        nist_id = int(l.strip().split()[1])
        entry_num += 1
        while nist_id > entry_num:
            NoDB.store("")
            entry_num += 1
        continue
    if l.startswith("Num peaks: "):
        entry = cStringIO.StringIO()
        num_peaks = int(l.strip().split()[2])
        for i in range(num_peaks):
            vals = f.readline().strip().split()
            # NIST seems to have shifted to a 10000 based normalization...
            # (but not after Lib2Nist treatment, so, nevermind).
            # print >> entry, "%s,%d" % (vals[0],int(vals[1])/10)
            print >> entry, "%s,%s" % (vals[0], vals[1])
        NoDB.store(entry.getvalue())
        entry.close()
f.close()
NoDB.close()
