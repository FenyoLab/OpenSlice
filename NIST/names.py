import sys

f = open(sys.argv[1])
curname = ""
curid = -1
for l in f:
    if l.startswith("Name: "):
        curname = l.strip()[6:]
    if l.startswith("ID: "):
        curid = int(l.strip()[4:])
        print "%d\t%s" % (curid, curname)
f.close()
