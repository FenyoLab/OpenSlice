import NoDB
import cStringIO

msps = ["cptac2_human_hcd_selected.MSP"]
msp_counter = 0
baseline = 0
f = open(msps[msp_counter])
NoDB.create("nist")
entry_num = 0
nist_id = 0
while True:
    l = f.readline()
    if not l:
        f.close()
        msp_counter += 1
        if msp_counter == len(msps):
            break
        f = open(msps[msp_counter])
        baseline = entry_num
        print msps[msp_counter], baseline
        l = f.readline()
    if l.startswith("ID: "):
        nist_id = int(l.strip().split()[1]) + baseline
        entry_num += 1
        while nist_id > entry_num:
            NoDB.store("")
            entry_num += 1
        continue
    if l.startswith("Num peaks: "):
        entry = cStringIO.StringIO()
        # Discovered some 'nominal mass' spectra in the MSPs, must guard from this 
        # and change decoding when this happens...
        #     if ";" in first_line:
        while True:
            a_line = f.readline().strip()
            if a_line == "":  # Assumes an empty entry even in the last line of the file...
                NoDB.store(entry.getvalue())
                entry.close()
                break
            if ";" in a_line:  # Assumes there is never a ';' in peak commentary...
                print a_line
                vals = a_line.split(";")
                if vals[-1] != '':
                    print "last entry of nominal mass MSP entry should have been empty..."
                    import sys

                    sys.exit(-1)
                for pair in vals[:-1]:
                    pair_vals = pair.split()
                    print >> entry, "%s,%s" % (pair_vals[0], pair_vals[1])
            else:  # This assumes that in nominal mass entries there is a trailing ';' even on a single-pair line.
                vals = a_line.split()
                print >> entry, "%s,%s" % (vals[0], vals[1])
NoDB.close()
