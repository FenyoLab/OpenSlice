import os

fnames = [f for f in os.listdir(".") if f.endswith(".mzdb")]
sample_names = [fname[:-5] for fname in fnames]
completed = [fname[:-14] for fname in os.listdir(".") if fname.endswith("_completed.txt")]
todo_samples = set(sample_names).difference(completed)


# This level of indirection exists only because of SQLite's
# weird behaviour -- there is some weird memory leak or something
# so I defensively launch a process per file. Sad, but necessary.
def zeno_launch(sample_name):
    os.system("python zeno.py %s" % sample_name)


from multiprocessing import Pool

if __name__ == "__main__":
    p = Pool(8)
    p.map(zeno_launch, todo_samples)
