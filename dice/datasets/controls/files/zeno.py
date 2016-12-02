import time
import cStringIO
import NoDB
import No_mz
from mz.API import mzdb
import sys
import os

sample_name = sys.argv[1]


def block_pass(min_nmz, max_nmz):
    # does not include max_smz (to avoid duplicate runs)...
    pass_start = time.time()

    xics = {}
    # nmz values will range from 1 to N-1 where N is supplied by No_mz.N... (see No_mz for details).
    for nmz in range(min_nmz, max_nmz):
        xics[nmz] = []

    t_counter = 0
    for t in times:
        if t_counter and not t_counter % 1000:
            print t_counter
        t_counter += 1
        xic_accs = {}
        for nmz in range(min_nmz, max_nmz):
            xic_accs[nmz] = 0.0
        # Can tighten this to only include relevant mass range, but since there is overlap between slices, this can
        # get a bit tricky. For the time being, I prefer to be somewhat wasteful but more verbose and clear...
        for (mz, y) in [(mz, y) for (mz, y, _s, z) in raw.cscan(raw.scanForTime(t)) if min_mz <= mz <= max_mz]:
            offsets = No_mz.numeros(mz)
            for offset in offsets:
                if min_nmz <= offset < max_nmz:
                    xic_accs[offset] += y
        for nmz in range(min_nmz, max_nmz):
            # XICs will always be rounded to the nearest integer...
            xics[nmz].append(int(xic_accs[nmz] + 0.5))

    for nmz in range(min_nmz, max_nmz):
        f = cStringIO.StringIO()
        put_minus = True
        counter = 1
        for v in xics[nmz]:
            if v == 0.0:
                put_minus = True
            else:
                if put_minus:
                    print >> f, -1 * counter
                    put_minus = False
                print >> f, v
            counter += 1
        entry_num = NoDB.store(f.getvalue())
        assert entry_num == nmz
        f.close()
    pass_stop = time.time()
    print "%d ICs generated in %.2f seconds" % (max_nmz - min_nmz, pass_stop - pass_start)


total_start = time.time()
if not os.path.exists(sample_name):
    os.mkdir(sample_name)

raw = mzdb.File(sample_name + ".mzdb")
NoDB.create(sample_name + "/xics")

ms1_times = dict([(t, offset + 1) for (offset, t) in enumerate([t for (t, p) in raw.scan_list() if p == 0.0])])
f = open(sample_name + "/xic_times.csv", 'w')
times = ms1_times.keys()
times.sort()
for t in times:
    print >> f, t
f.close()

print "ms1: %d" % len(ms1_times)

(min_mz, max_mz) = No_mz.domain
# Hark -- a magic number (25000) !?! 
block_size = 25000
boundaries = []
bottom = 1
while bottom < No_mz.N:
    boundaries.append((bottom, min(No_mz.N, bottom + block_size)))
    bottom += block_size

for (bottom, top) in boundaries:
    print "%d-%d (of %d)" % (bottom, top, No_mz.N - 1)
    block_pass(bottom, top)

NoDB.close()
raw.close()
total_stop = time.time()
print "Rawfile extracted in %.2f seconds" % (total_stop - total_start)

output = open(sample_name + "_completed.txt", 'w')
print >> output, "Rawfile extracted in %.2f seconds" % (total_stop - total_start)
output.close()
