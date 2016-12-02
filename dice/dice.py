import NoDB
from bottle import route, get, static_file, abort, run

#
# Will implement log scale...
#

#
# Initalization...
#
datasets = set()
f = open("datasets/published.txt")
for l in f:
    datasets.add(l.strip())
f.close()
# To avoid a project named 'static' from being clobbered by the
# webserver override on the path "/dice/static/"...
assert not "static" in datasets
assert not "nist" in datasets


#
#  Common Static Files... (dataset variable ignored)
#
@route('/<dataset>/ms2_full_viewer.html')
def ms2_full_viewer(dataset):
    return static_file("ms2_full_viewer.html", "static")


@route('/<dataset>/ms2_viewer.html')
def ms2_viewer(dataset):
    return static_file("ms2_viewer.html", "static")


@route('/<dataset>/nist_viewer.html')
def ms2_nist_viewer(dataset):
    return static_file("nist_viewer.html", "static")


@get('/<dataset>/js/<filename:re:.*\.js>')
def dataset_javascripts(dataset, filename):
    return static_file(filename, root='static/js')


@get('/<dataset>/css/<filename:re:.*\.css>')
def dataset_stylesheets(dataset, filename):
    return static_file(filename, root='static/css')


@get('/<dataset>/img/<filename:re:.*\.(jpg|png|gif|ico)>')
def dataset_images(dataset, filename):
    return static_file(filename, root='static/img')


@get('/<dataset>/font/<filename:re:.*\.(eot|ttf|woff|svg)>')
def dataset_fonts(dataset, filename):
    return static_file(filename, root='static/fonts')


@route('/')
def index():
    return static_file("index.html", "static")


#
#  Supporting datasystem
#
@route('/nist/<no:int>')
def nist(no):
    return NoDB.lookup("datasets/nist", no)


#
#  Dataset dependent -- first check that the path
#                       contains a published dataset
#
@route('/<dataset>/files')
def files(dataset):
    if dataset in datasets:
        return static_file("published.txt", "datasets/" + dataset)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


@route('/<dataset>/published')
def published(dataset):
    if dataset in datasets:
        return static_file("published.json", "datasets/" + dataset)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


@route('/<dataset>/xic_times')
def pic_times(dataset):
    if dataset in datasets:
        return static_file("xic_times.json", "datasets/" + dataset)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


import os
import os.path
import json


@route('/<dataset>/pic/<no:int>')
def pic(dataset, no):
    if dataset in datasets:
        if os.path.exists("datasets/" + dataset + "/pics/%d.json" % no):
            return static_file("%d.json" % no, "datasets/" + dataset + "/pics")
        else:
            the_pic = {}
            f = open("datasets/%s/published.txt" % dataset)
            fnames = [l.strip() for l in f]
            f.close()
            print fnames
            for fn in fnames:
                the_pic[fn] = map(int, NoDB.lookup("datasets/" + dataset + "/files/" + fn + "/xics", no).split())
                print the_pic[fn]
            out = open("datasets/" + dataset + "/pics/%d.json" % no, 'w')
            print >> out, json.dumps(the_pic)
            out.close()
            return static_file("%d.json" % no, "datasets/" + dataset + "/pics")
    else:
        abort(404, "Invalid dataset (%s)." % datasets)


@route('/<dataset>/precursors/<no:int>')
def precursors(dataset, no):
    if dataset in datasets:
        return NoDB.lookup("datasets/" + dataset + "/precursors", no)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


@route('/<dataset>/files/<fname>/xics/<no:int>')
def xic(dataset, fname, no):
    if dataset in datasets:
        return NoDB.lookup("datasets/" + dataset + "/files/" + fname + "/xics", no)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


@route('/<dataset>/files/<fname>/scans/<no:int>')
def scan(dataset, fname, no):
    if dataset in datasets:
        return NoDB.lookup("datasets/" + dataset + "/files/" + fname + "/scans", no)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


@route('/<dataset>/files/<fname>/xic_times.csv')
def xic_times(dataset, fname):
    if dataset in datasets:
        return static_file("xic_times.csv", "datasets/" + dataset + "/files/" + fname)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


#
#  Static, but dataset specific...
#
import os


@route('/<dataset>/<page>.html')
def dataset_specific_page(dataset, page):
    if dataset in datasets:
        if page + ".html" in os.listdir("datasets/" + dataset + "/views"):
            return static_file(page + ".html", "datasets/" + dataset + "/views")
        else:
            abort(404, "page %s.html not found in dataset %s." % (page, dataset))
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


@route('/<dataset>/')
def dataset_index(dataset):
    if dataset in datasets:
        if "index.html" in os.listdir("datasets/" + dataset + "/views"):
            return static_file("index.html", "datasets/" + dataset + "/views")
        else:
            abort(404, "page index.html not found in dataset %s." % dataset)
    else:
        abort(404, "Invalid dataset (%s)." % dataset)


#
# Default Static Routes
#
@get('/js/<filename:re:.*\.js>')
def javascripts(filename):
    return static_file(filename, root='static/js')


@get('/css/<filename:re:.*\.css>')
def stylesheets(filename):
    return static_file(filename, root='static/css')


@get('/img/<filename:re:.*\.(jpg|png|gif|ico)>')
def images(filename):
    return static_file(filename, root='static/img')


@get('/font/<filename:re:.*\.(eot|ttf|woff|svg)>')
def fonts(filename):
    return static_file(filename, root='static/fonts')


@get('/ptpatlas/<protein>')
def ptpatlas(protein):
    return static_file(protein, root='../ptpatlas')


@get('/proteinset/<pset>')
def proteinset(pset):
    return static_file(pset, root='../proteinset')


if __name__ == "__main__":
    run(host='0.0.0.0', port=80, debug=True)
