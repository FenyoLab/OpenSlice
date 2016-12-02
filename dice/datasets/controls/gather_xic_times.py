xic_times = {}
f = open("published.txt")
for l in f:
    fname = l.strip()
    tf = open("files/" + fname + "/xic_times.csv")
    times = []
    for lf in tf:
        times.append(float(lf.strip()))
    tf.close()
    xic_times[fname] = times
f.close()

import json

out = open("xic_times.json", 'w')
print >> out, json.dumps(xic_times)
out.close()
