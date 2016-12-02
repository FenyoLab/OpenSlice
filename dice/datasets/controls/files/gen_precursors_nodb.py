import cStringIO
import NoDB
import No_mz
import sys
import os

todo_samples = [fname[:-14] for fname in os.listdir(".") if fname.endswith("_completed.txt")]
scans = {}
for sample in todo_samples:
    f = open(sample + "/scans.tab")
    for l in f:
        (stime, truval, z, stype, snum, nist_id, peptide, nist_score) = l.strip("\n").split("\t")
        if truval == "0.0000":
            continue
        # Hopefully this causes no warpage in mz-value...
        mskeys = No_mz.numeros(float(truval))
        for mskey in mskeys:
            if not mskey in scans:
                scans[mskey] = []
            scans[mskey].append((sample, stime, z, stype, truval, snum, nist_id, peptide, nist_score))
    f.close()
print len(scans)

NoDB.create("precursors")
for k in range(1, No_mz.N):
    entry = ""
    if k in scans:
        f = cStringIO.StringIO()
        the_scans = scans[k]
        the_scans.sort()
        for (sample, stime, z, stype, truval, snum, nist_id, peptide, nist_score) in scans[k]:
            print >> f, "\t".join([sample, stime, z, stype, truval, snum, nist_id, peptide, nist_score])
        entry = f.getvalue()
        f.close()
    entry_num = NoDB.store(entry)
    assert entry_num == k
NoDB.close()
