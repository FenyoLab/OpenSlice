import os

to_publish = [fname[:-14] for fname in os.listdir("files") if fname.endswith("_completed.txt")]

f = open("published.txt", 'w');
for sample in to_publish:
    print >> f, sample
f.close()
