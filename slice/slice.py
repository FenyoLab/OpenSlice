import sys

sys.path.append(".")

import os
import json

# DO THIS BEFORE importing mz!!!
import tempfile

os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib

import mz
import mz.tools
from mz.API import mzdb

from bottle import route, run, response, get, static_file

# import matplotlib
# import matplotlib.cbook
# from mz import tools
from mako.template import Template

mimetypes = {"html": "text/html", "csv": "text/html", "json": "application/json", "txt": "text/plain",
             "png": "image/png", "svg": "image/svg+xml", "mgf": "text/plain", "dta": "text/plain",
             "pdf": "application/pdf"}


# TODO: rewrite this using regexps later...
def untype(aString):
    for m in mimetypes:
        if aString.endswith("." + m):
            typelen = len(m) + 1
            return (aString[:-typelen], m)
    return (aString, "html")


def to_last_slash(aString):
    bString = aString[::-1]
    return aString[:-1 * bString.find("/")]


# Short Float...
# TODO: I know there is a better way... Fix this later...
def shf(aFloat, k=3):
    format = "%%.%df" % k
    return format % aFloat


def svg_match_request_url(filename, scan_time, sequence):
    return baseurl + "/files/%s/scans/%s/match/%s.svg" % (filename, shf(scan_time), sequence)


def png_match_request_url(dataset, filename, scan_time, sequence):
    return "scans/%s/match/%s.png" % (shf(scan_time), sequence)


def png_match_request_url_from_match(dataset, filename, scan_time, sequence):
    return "../../../scans/%s/match/%s.png" % (shf(scan_time), sequence)


def ric_request_url(dataset, filename, time_start, time_stop, mz_start, mz_stop):
    return "../xic/%s-%s/%s-%s" % (shf(time_start), shf(time_stop), shf(mz_start), shf(mz_stop))


def ric_request_url_from_overlay(dataset, filename, time_start, time_stop, mz_start, mz_stop):
    return "../../../xic/%s-%s/%s-%s" % (shf(time_start), shf(time_stop), shf(mz_start), shf(mz_stop))


def scan_request_url(dataset, filename, scan_time):
    return "%s" % shf(scan_time)


def scan_request_url_from_overlay(dataset, filename, scan_time):
    return "../../%s" % shf(scan_time)


def generic_scan_request_url(dataset, filename):
    return "../../scans/"


def pic_template(dataset, pic_report, mz, tstart, tstop):
    pic_template = Template(filename='../dice/datasets/' + dataset + '/views/pic.html')
    return pic_template.render(pic_report=pic_report, tmz=mz, tstart=tstart, tstop=tstop)


def ric_json_template(dataset, ric_data, mz_center, ctype, filename):
    data = "t,intensity,msms_intensity\\n"
    for (t, s, msms_s) in ric_data:
        if msms_s:
            data += "%.3f,%.0f,%.0f\\n" % (t, s, msms_s)
        else:
            data += "%.3f,%.0f,\\n" % (t, s)
    rurl_link = generic_scan_request_url(dataset, filename)
    ric_template = Template(filename='views/ric.html')

    mz_center_str = "%.4f" % mz_center
    return ric_template.render(data=data, ctype=ctype, mz_center=mz_center_str, rurl_link=rurl_link)


def bpc_json_template(dataset, bpc_data, ctype, filename):
    data = "t,intensity,msms_intensity\\n"
    for (t, s, msms_s) in bpc_data:
        if msms_s:
            data += "%.3f,%.0f,%.0f\\n" % (t, s, msms_s)
        else:
            data += "%.3f,%.0f,\\n" % (t, s)
    rurl_link = filename + "/scans/"
    bpc_template = Template(filename='views/bpc.html')
    return bpc_template.render(data=data, ctype=ctype, rurl_link=rurl_link)


def ric_template(dataset, ric_data, mz_center, ctype, filename):
    if ctype == "html":
        return ric_json_template(dataset, ric_data, mz_center, ctype, filename)

    from cStringIO import StringIO

    xic_time = []
    xic_int = []
    scan_dot = None  # I don't currently distinguish the "source" scan...
    bin_times = []
    bin_ints = []
    for (t, i, ii) in ric_data:
        if ii:
            bin_times.append(t)
            bin_ints.append(ii)
        else:
            xic_time.append(t)
            xic_int.append(i)
    tmp = StringIO()
    if ctype == "png":
        fmt = "PNG"
    if ctype == "svg":
        fmt = "SVG"
    if ctype == "pdf":
        fmt = "PDF"

    make_xic_im(tmp, mz_center, xic_time, xic_int, scan_dot, bin_times, bin_ints, fmt)

    text = tmp.getvalue()

    return text


def _make_xic(fig, xic_center, xic_time, xic_int, scan_dot, bin_times, bin_ints):
    '''Plots a XIC, with labels for the primary MSMS scan and any neighboring scans
    in the same time and mz area.

    This is an internal method, it expects a matplotlib Figure instance'''

    from matplotlib.ticker import ScalarFormatter

    fig.clear()
    fig.set_facecolor('w')

    axes = fig.add_axes([0.125, 0.1, 0.775, 0.8])

    axes.plot(xic_time, xic_int, '-g', linewidth=2, markeredgecolor='g',
              markerfacecolor='g', markersize=2)

    if len(bin_times) > 0:
        axes.plot(bin_times, bin_ints, 'bo', markersize=4)

    # plot the scan dot last, so it stays on top of other markers
    if scan_dot is not None:
        axes.plot([scan_dot[0]], [scan_dot[1]], 'b^', markersize=10)

    axes.set_xlabel('Time (min)')
    axes.set_ylabel('Abundance')

    axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

    return (zip(xic_time, xic_int) + ([scan_dot] if scan_dot else []),)


def _make_bpc(fig, xic_time, xic_int, bin_times, bin_ints):
    '''Plots a BPC, with labels for the primary MSMS scan.

    This is an internal method, it expects a matplotlib Figure instance'''

    from matplotlib.ticker import ScalarFormatter

    fig.clear()
    fig.set_facecolor('w')

    axes = fig.add_axes([0.125, 0.1, 0.775, 0.8])

    axes.plot(xic_time, xic_int, '-g', linewidth=2, markeredgecolor='g',
              markerfacecolor='g', markersize=2)

    if len(bin_times) > 0:
        axes.plot(bin_times, bin_ints, 'bo', markersize=4)

    axes.set_title('Base Peak Chromatogram')
    axes.set_xlabel('Time (min)')
    axes.set_ylabel('Abundance')

    axes.xaxis.set_major_formatter(ScalarFormatter(useOffset=False, useMathText=True))

    return (zip(xic_time, xic_int),)


def make_xic_im(save_file, xic_center, xic_time, xic_int, scan_dot, bin_times, bin_ints, fmt='PNG'):
    '''Plots an XIC, with labels for the primary MSMS scan and any neighboring scans
    in the same time and mz area'''
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig = Figure()
    _make_xic(fig, xic_center, xic_time, xic_int, scan_dot, bin_times, bin_ints)
    FigureCanvasAgg(fig).print_figure(save_file, dpi=80, format=fmt)


def make_bpc_im(save_file, xic_time, xic_int, bin_times, bin_ints, fmt='PNG'):
    '''Plots a BPC, with labels for the primary MSMS scans'''
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    fig = Figure()
    _make_bpc(fig, xic_time, xic_int, bin_times, bin_ints)
    FigureCanvasAgg(fig).print_figure(save_file, dpi=80, format=fmt)


def bpc_template(dataset, bpc_data, ctype, filename):
    if ctype == "html":
        return bpc_json_template(dataset, bpc_data, ctype, filename)

    from cStringIO import StringIO

    xic_time = []
    xic_int = []
    bin_times = []
    bin_ints = []
    for (t, i, ii) in bpc_data:
        if ii:
            bin_times.append(t)
            bin_ints.append(ii)
        else:
            xic_time.append(t)
            xic_int.append(i)
    tmp = StringIO()
    if ctype == "png":
        fmt = "PNG"
    if ctype == "svg":
        fmt = "SVG"
    if ctype == "pdf":
        fmt = "PDF"

    make_bpc_im(tmp, xic_time, xic_int, bin_times, bin_ints, fmt)

    text = tmp.getvalue()

    return text


def scan_json_template(dataset, mzScan, scan_time, ctype, filename, last_ms1, prev_time, ric_mass, next_time, next_ms1,
                       charge, scan_num):
    if ric_mass == 0.0:
        disabled = 'disabled="disabled"'
        ico = "ms1ico.png"
        title = "time: %.3f min (scan# %d)" % (float(scan_time), scan_num)
    else:
        chargestr = ""
        if charge:
            chargestr = "(%d+)" % charge
        scanstr = ""
        if scan_num:
            scanstr = " (scan# %d)" % scan_num

        disabled = ""
        ico = "ms2ico.png"
        title = "m/z=%.3f%s &nbsp;&nbsp; time: %.2fmin%s" % (ric_mass, chargestr, float(scan_time), scanstr)

    data = "mz,intensity\\n"
    zs = ""
    for vals in mzScan:
        if len(vals) == 2:
            data += "%f,%f\\n" % (vals[0], vals[1])
        else:
            data += "%f,%f\\n" % (vals[0], vals[1])
            zs += "%d\\n" % vals[2]

    scan_template = Template(filename='views/scan.html')
    prev_ms1_link = scan_request_url(dataset, filename, last_ms1)
    prev_scan_link = scan_request_url(dataset, filename, prev_time)
    ric_link = ric_request_url(dataset, filename, float(scan_time) - 3, float(scan_time) + 3, ric_mass - 0.05,
                               ric_mass + 0.05)
    next_scan_link = scan_request_url(dataset, filename, next_time)
    next_ms1_link = scan_request_url(dataset, filename, next_ms1)
    ric_last_slash = to_last_slash(ric_link)
    project_link = "../../../"
    text = scan_template.render(ico=ico, title=title, disabled=disabled, data=data, zs=zs,
                                ric_last_slash=ric_last_slash, ctype=ctype, prev_ms1_link=prev_ms1_link,
                                prev_scan_link=prev_scan_link, ric_link=ric_link, next_scan_link=next_scan_link,
                                next_ms1_link=next_ms1_link, project_link=project_link,
                                scan_time="%.2f" % float(scan_time))

    return text


def scan_template(dataset, mzScan, scan_time, ctype, filename, last_ms1, prev_time, ric_mass, next_time, next_ms1,
                  charge, scan_num):
    if ctype == "html":
        return scan_json_template(dataset, mzScan, scan_time, ctype, filename, last_ms1, prev_time, ric_mass, next_time,
                                  next_ms1, charge, scan_num)

    if ctype == "png":
        fmt = "PNG"
    if ctype == "svg":
        fmt = "SVG"
    if ctype == "pdf":
        fmt = "PDF"

    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from cStringIO import StringIO

    precursor_ions = []
    for k in range(1, charge + 1):
        precursor_ions.append("MH%s" % "".join(["+"] * k))

    fig = Figure()
    tmp = StringIO()

    if ric_mass == 0:
        mzScan = [(x, y) for (x, y, z) in mzScan]
        charge = 0
        precursor_mz = 0
        mz.tools.mz_image._make_ms2(fig, mzScan, "c", None,
                                    title="time=%.2fmin  (scan# %d)" % (float(scan_time), scan_num), labels=None,
                                    ion_list=(['b', 'y'] + precursor_ions), filter_labels=False, charge=charge)
    else:
        chargestr = ""
        if charge:
            chargestr = "(%d+)" % charge
        scanstr = ""
        if scan_num:
            scanstr = " (scan# %d)" % scan_num
        print precursor_ions
        mz.tools.mz_image._make_ms2(fig, mzScan, "c", None, labels=None, ion_list=(['b', 'y'] + precursor_ions),
                                    title="m/z=%.3f%s  time=%.2fmin%s" % (
                                    float(ric_mass), chargestr, float(scan_time), scanstr), filter_labels=False,
                                    charge=charge)
    FigureCanvasAgg(fig).print_figure(tmp, format=fmt, dpi=80)
    text = tmp.getvalue()

    return text


def match_template(dataset, scan_time, mzScan, sequence, ctype, filename, last_ms1, prev_time, ric_mass, next_time,
                   next_ms1, scan_type="cid", charge=0, scan=0, precursor=0.0):
    text = ""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from cStringIO import StringIO
    from mz.functions import generate_labels

    fig = Figure()
    tmp = StringIO()
    chargestr = ""
    if charge:
        chargestr = "(%d+)" % charge
    scanstr = ""
    if scan:
        scanstr = " (scan# %d)" % scan

    if scan_type == "etd":
        if charge > 2:
            ion_list = ['c', 'z+1', 'c++', 'z+1++', 'MH']
        if charge > 0:
            ion_list = ['c', 'z+1', 'MH']
        if charge == 0:
            ion_list = ['c', 'z+1']
    else:
        if charge > 2:
            ion_list = ['b', 'b++', 'y', 'y++', 'MH']
        if charge > 0:
            ion_list = ['b', 'y', 'MH']
        if charge == 0:
            ion_list = ['b', 'y']

    filtered = [(an_mz, an_int) for (an_mz, an_int) in mzScan]
    filtered = sorted(filtered, key=lambda mi: mi[1], reverse=True)[:50]
    labels = generate_labels(filtered, sequence, ions=ion_list, tolerance=0.6, charge=charge)

    new_labels = []
    for (x, lab) in labels:
        new_labels.append((x, lab.split()[0]))
    (labeled, ignore) = mz.tools.mz_image._make_ms2(fig, mzScan, "c", sequence, title="m/z=%.3f%s  time=%.2fmin%s" % (
    float(ric_mass), chargestr, float(scan_time), scanstr), labels=new_labels, ion_list=ion_list, charge=charge,
                                                    filter_labels=False)
    if ctype in ["pdf", "png", "svg"]:
        img_format = ctype
        FigureCanvasAgg(fig).print_figure(tmp, format=img_format[-3:].upper(), dpi=80)
    if ctype in ["png", "svg", "pdf"]:
        text = tmp.getvalue()
    if ctype == "svg":
        # pass
        text = text.replace('viewBox="0 0 576 432"', '')
        text = text.replace('<g id="figure_1">',
                            '<script xlink:href="../../SVGPan.js"/><g id="viewport" transform="scale(1.12)">')
    if ctype == "html":
        if ric_mass == 0.0:
            disabled = 'disabled="disabled"'
        else:
            disabled = ""

        match_template = Template(filename='views/match.html')
        scan_time_str = "%.2f" % float(scan_time)
        pnglink = png_match_request_url_from_match(dataset, filename, float(scan_time), sequence)
        prev_ms1_link = scan_request_url(dataset, filename, last_ms1)
        prev_scan_link = scan_request_url(dataset, filename, prev_time)
        ric_link = ric_request_url(dataset, filename, float(scan_time) - 3, float(scan_time) + 3, ric_mass - 0.05,
                                   ric_mass + 0.05)
        next_scan_link = scan_request_url(dataset, filename, next_time)
        next_ms1_link = scan_request_url(dataset, filename, next_ms1)

        text = match_template.render(scan_time=scan_time_str, disabled=disabled, prev_ms1_link=prev_ms1_link,
                                     prev_scan_link=prev_scan_link, ric_link=ric_link, next_scan_link=next_scan_link,
                                     next_ms1_link=next_ms1_link, pnglink=pnglink)

    if ctype == ".txt":
        text = `labeled`

    return text


mod_lookup = {"c": 57.021464, "m": 15.994915, "p": 15.994915, "s": 79.966331, "t": 79.966331, "y": 79.966331}


def overlay_template(dataset, scan_time, mzScan, sequence, ctype, filename, last_ms1, prev_time, ric_mass, next_time,
                     next_ms1, scan_type="cid", charge=0, scan=0, precursor=0.0, format="minimal"):
    text = ""

    chargestr = ""
    if charge:
        chargestr = "(%d+)" % charge
    scanstr = ""
    if scan:
        scanstr = " (scan# %d)" % scan

    if ric_mass == 0.0:
        disabled = 'disabled="disabled"'
    else:
        disabled = ""

    if format == "minimal":
        overlay_template = Template(filename='views/overlay.html')
    else:
        overlay_template = Template(filename='views/full_overlay.html')

    scan_time_str = "%.2f" % float(scan_time)
    prev_ms1_link = scan_request_url_from_overlay(dataset, filename, last_ms1)
    prev_scan_link = scan_request_url_from_overlay(dataset, filename, prev_time)
    ric_link = ric_request_url_from_overlay(dataset, filename, float(scan_time) - 3, float(scan_time) + 3,
                                            ric_mass - 0.05, ric_mass + 0.05)
    next_scan_link = scan_request_url_from_overlay(dataset, filename, next_time)
    next_ms1_link = scan_request_url_from_overlay(dataset, filename, next_ms1)

    scan_struct = {}
    variable_mods = []
    for (offset, letter) in enumerate(sequence):
        if letter.islower():
            variable_mods.append({"index": offset + 1, "modMass": mod_lookup[letter], "aminoAcid": letter.upper()})
    scan_struct["sequence"] = sequence.upper()
    scan_struct["scanNum"] = scan
    scan_struct["charge"] = charge
    scan_struct["precursorMz"] = precursor
    scan_struct["variableMods"] = variable_mods
    scan_struct["fileName"] = filename
    scan_struct["peaks"] = [[t[0], t[1]] for t in mzScan]
    scan_struct["format"] = format

    text = overlay_template.render(scan_time=scan_time_str, disabled=disabled, prev_ms1_link=prev_ms1_link,
                                   prev_scan_link=prev_scan_link, ric_link=ric_link, next_scan_link=next_scan_link,
                                   next_ms1_link=next_ms1_link, lorikeet_json=json.dumps(scan_struct, indent=2),
                                   overlayurl="/files/" + filename + "/scans/" + shf(float(scan_time)) + "/overlay/",
                                   sequence=sequence)

    return text


def files(dataset):
    import glob
    mzdbfiles = glob.glob("../dice/datasets/%s/files/*.mzdb" % dataset)
    files = mzdbfiles
    files = map(os.path.basename, files)
    files = [os.path.splitext(s)[0] for s in files]
    files.sort()
    return files


cached = []
cached_dataset = ""


def check_cache(dataset, target):
    global cached_dataset
    global cached
    if cached_dataset != dataset:
        cached_dataset = dataset
        cached = []
        if os.path.exists("../dice/datasets/" + dataset + "/cache/cached.txt"):
            f = open("../dice/datasets/" + dataset + "/cache/cached.txt", 'U')
            for l in f:
                print l
                cached.append(float(l.strip()))
            f.close()
            cached.sort()
        print cached
    best_diff = 1000000.0
    last_delta = best_diff
    hit = None
    for entry in cached:
        delta = abs(entry - float(target))
        if delta > last_delta:
            break
        if delta < best_diff:
            best_diff = delta
            hit = entry
        last_delta = delta
    if hit and abs(hit - float(target)) / float(target) < 5.0 / 1000000.0:
        print "**** Cache entry found!!! ****"
        f = open("../dice/datasets/" + dataset + "/cache/%.4f" % hit)
        output = f.read()
        f.close()
        return output
    else:
        ric_data = {}
        start_mz_float = float(target) * (1.0 - 5E-6)
        stop_mz_float = float(target) * (1.0 + 5E-6)
        for f in files(dataset):
            file = confirm_and_extract(dataset, f)
            ric_data[f] = file.xic(0.0, 1000.0, start_mz_float, stop_mz_float)
        file_list = ric_data.keys()
        file_list.sort()
        ric_report = {"ric_data": ric_data, "file_list": file_list, "memset": None, "metabolite": None,
                      "mz": "%.4f" % float(target), "rt": None}
        out = json.dumps(ric_report)
        f = open("../dice/datasets/" + dataset + "/cache/cached.txt", 'a')
        print >> f, "%.4f" % float(target)
        f.close()
        f = open("../dice/datasets/" + dataset + "/cache/%.4f" % float(target), 'w')
        print >> f, out
        f.close()
        cached.append(float(target))
        cached.sort()
        return out


def confirm_and_extract(dataset, filename):
    if os.path.exists("../dice/datasets/" + dataset + "/files/" + filename + ".mzdb"):
        print "got it!"
        retval = mzdb.File("../dice/datasets/" + dataset + "/files/" + filename + ".mzdb")
        print retval
        return retval
    else:
        print "don't got it!!!"
        return None


def do_scan_number_request(filename, scan, ctype):
    scan = int(scan)
    file = confirm_and_extract(filename)
    the_time = file.scan_list()[scan - 1][0]
    return do_scan_request(filename, the_time, ctype)


def do_scan_request(dataset, filename, scan_time, ctype):
    text = 0
    file = confirm_and_extract(dataset, filename)
    if file.file_type == "raw":
        the_filters = file.filters()
        (first, last) = file.scan_range()
        scan_num = file.scanForTime(float(scan_time))
        precursor = 0.0
        charge = 0
        if the_filters[scan_num - 1][1].find("FTMS") > -1:
            if ctype != ".txt":
                mzScan = [(mz, s, z) for (mz, s, n, z) in file.lscan(scan_num)]
            else:
                mzScan = [(mz, s) for (mz, s, n, z) in file.lscan(scan_num)]
        else:
            mzScan = file.scan(float(scan_time))
        try:
            precursor = 0.0
            charge = 0
            (precursor, charge) = file.scanPrecursor(scan_num)
        except:
            pass

        if ctype == "txt":
            text = ""
            for vals in mzScan:
                text += "%.6f\t%.0f\n" % (vals[0], vals[1])
        if ctype == "dta":
            text = ""
            if charge:
                text += "%.4f %d\n" % (charge * precursor - (charge - 1) * 1.007282675, charge)
            for vals in mzScan:
                if vals[1] >= 1.0:
                    text += "%.4f %.0f\n" % (vals[0], vals[1])
        if ctype == "mgf":
            text = ""
            text += "MASS=Monoisotopic\n"
            text += "SEARCH=MIS\n"
            text += "BEGIN IONS\n"
            text += "TITLE=%s\n" % scan_request_url(filename, float(scan_time))
            if precursor > 0.0:
                text += "PEPMASS=%.4f\n" % precursor
            if charge > 0:
                text += "CHARGE=%d+\n" % charge
            for vals in mzScan:
                if vals[1] >= 1.0:
                    text += "%.4f\t%.0f\n" % (vals[0], vals[1])
            text += "END IONS\n"

        if ctype not in ["txt", "mgf", "dta"]:
            if scan_num == last:
                next_time = scan_time
            else:
                next_time = file.timeForScan(scan_num + 1)
                if scan_num == first:
                    prev_time = scan_time
                else:
                    prev_time = file.timeForScan(scan_num - 1)
            if the_filters[scan_num - 1][1].find("Full ms ") > -1:
                print "I think this is an MS1..."
                ric_mass = 0.0
            else:
                print "I think this is an MS2..."
                try:
                    ric_mass = file.scanPrecursor(scan_num)[0]
                except:
                    ric_mass = 0.0
            print "RIC MASS:", ric_mass
            next_ms1 = -1.0
            ms1_nscan = scan_num
            while ms1_nscan != last:
                ms1_nscan += 1
                if the_filters[ms1_nscan - 1][1].find("Full ms ") > -1:
                    next_ms1 = file.timeForScan(ms1_nscan)
                    break
            last_ms1 = -1.0
            ms1_lscan = scan_num
            while ms1_lscan != first:
                ms1_lscan -= 1
                if the_filters[ms1_lscan - 1][1].find("Full ms ") > -1:
                    last_ms1 = file.timeForScan(ms1_lscan)
                    break

            text = scan_template(dataset, mzScan, scan_time, ctype, filename, last_ms1, prev_time, ric_mass, next_time,
                                 next_ms1, charge, scan_num)
    return text


def do_header_request(dataset, filename, scan_time, ctype):
    file = confirm_and_extract(dataset, filename)
    attributes = file.extra_info(file.scanForTime(float(scan_time)))
    if ctype == "html":
        text = "<html><head></head><body>"
        text += "<table>\n"
        for (attr, val) in attributes.items():
            text += "<tr>"
            text += '<td align="right">'
            text += "%s:" % attr
            text += "</td>"
            text += "<td>"
            if type(val) == float:
                strval = "%.4f" % val
            else:
                strval = `val`
            text += "%s" % strval
            text += "</td>"
            text += "</tr>\n"
        text += "</table></body></html>"
    else:
        text = ""
        for (attr, val) in attributes.items():
            if type(val) == float:
                strval = "%.4f" % val
            else:
                strval = `val`
            text += "%s\t%s\n" % (attr, strval)
    return text


def do_filter_request(dataset, filename, scan_time, ctype):
    file = confirm_and_extract(dataset, filename)
    filter = file.filters()[file.scanForTime(float(scan_time)) - 1][1]
    if ctype == "html":
        text = "<html><head></head><body>"
        text += filter
        text += "</body></html>"
    else:
        text = filter
    return text


def do_mic_request(dataset, filename, start_time, stop_time, start_mz, stop_mz, ctype):
    text = ""
    file = confirm_and_extract(dataset, filename)
    start_time_float = float(start_time)
    stop_time_float = float(stop_time)
    start_mz_float = float(start_mz)
    stop_mz_float = float(stop_mz)
    ric_data = file.xic(start_time_float, stop_time_float, start_mz_float, stop_mz_float)
    msmsts = [t for (t, mz, sname, stype, smode) in
              file.scan_info(start_time_float, stop_time_float, start_mz_float - 0.01, stop_mz_float + 0.01) if
              stype != "MS1"]
    msmsts.sort()
    prev_t = start_time_float
    prev_s = 0.0
    aug_ric_data = []
    for (t, s) in ric_data:
        while len(msmsts) and msmsts[0] > prev_t and msmsts[0] < t:
            inter = msmsts[0]
            msmsts = msmsts[1:]
            prop = (inter - prev_t) / (t - prev_t)
            inter_s = prop * s + (1 - prop) * prev_s
            aug_ric_data.append((inter, inter_s, inter_s))
        aug_ric_data.append((t, s, None))
        prev_t = t
        prev_s = s
    t = stop_time_float
    s = 0.0
    for inter in msmsts:
        prop = (inter - prev_t) / (t - prev_t)
        inter_s = prop * s + (1 - prop) * prev_s
        aug_ric_data.append((inter, inter_s, inter_s))

    if ctype == "txt":
        text = ""
        for (t, s, msms_s) in aug_ric_data:
            text += "%.3f\t%.0f\n" % (t, s)
    if ctype == "csv":
        text = ""
        for (t, s, msms_s) in aug_ric_data:
            if msms_s:
                text += "%.3f,%.0f,%.0f\n" % (t, s, msms_s)
            else:
                text += "%.3f,%.0f,\n" % (t, s)
    if ctype == "json":
        text = "["
        for (t, s, msms_s) in aug_ric_data[:-1]:
            if msms_s:
                text += "[%.3f,%.0f,1],\n" % (t, s)
            else:
                text += "[%.3f,%.0f,0],\n" % (t, s)
        (t, s, msms_s) = aug_ric_data[-1]
        if msms_s:
            text += "[%.3f,%.0f,1]]\n" % (t, s)
        else:
            text += "[%.3f,%.0f,0]]\n" % (t, s)
    if ctype not in ["txt", "csv", "json"]:
        text = ric_template(dataset, aug_ric_data, (start_mz_float + stop_mz_float) / 2.0, ctype, filename)
    return text


def do_jpic_request(dataset, mz, rt):
    out = check_cache(dataset, mz)
    if out:
        return out
    ric_data = {}
    for f in files(dataset):
        file = confirm_and_extract(dataset, f)
        start_time_float = float(rt) - 0.5
        stop_time_float = float(rt) + 0.5
        start_mz_float = float(mz) * (1.0 - 5E-6)
        stop_mz_float = float(mz) * (1.0 + 5E-6)
        ric_data[f] = file.xic(0.0, 1000.0, start_mz_float, stop_mz_float)
    file_list = ric_data.keys()
    file_list.sort()
    ric_report = {"ric_data": ric_data, "file_list": file_list, "memset": None, "metabolite": None, "mz": mz, "rt": rt}
    return json.dumps(ric_report)


def do_pic_request(dataset, mz, rt):
    start_time_float = float(rt) - 0.5
    stop_time_float = float(rt) + 0.5
    start_mz_float = float(mz) * (1.0 - 5E-6)
    stop_mz_float = float(mz) * (1.0 + 5E-6)
    out = check_cache(dataset, mz)
    if out:
        return pic_template(dataset, out, mz, start_time_float, stop_time_float)
    text = ""
    ric_data = {}
    for f in files(dataset):
        file = confirm_and_extract(dataset, f)
        ric_data[f] = file.xic(0.0, 1000.0, start_mz_float, stop_mz_float)
    file_list = ric_data.keys()
    file_list.sort()
    ric_report = {"ric_data": ric_data, "file_list": file_list, "memset": None, "metabolite": None, "mz": mz, "rt": rt}
    text = pic_template(dataset, json.dumps(ric_report), mz, start_time_float, stop_time_float)
    return text


def do_bpc_request(dataset, filename, ctype):
    text = ""
    file = confirm_and_extract(dataset, filename)
    (start_time_float, stop_time_float) = file.time_range()
    bpc_data = file.bpc()
    # Twilight
    msmsts = [t for (t, mz, sname, stype, smode) in file.scan_info(start_time_float, stop_time_float, 0.0, 1000000.0) if
              stype != "MS1"]
    msmsts.sort()
    prev_t = start_time_float
    prev_s = 0.0
    aug_bpc_data = []
    for (t, s) in bpc_data:
        while len(msmsts) and msmsts[0] > prev_t and msmsts[0] < t:
            inter = msmsts[0]
            msmsts = msmsts[1:]
            prop = (inter - prev_t) / (t - prev_t)
            inter_s = prop * s + (1 - prop) * prev_s
            aug_bpc_data.append((inter, inter_s, inter_s))
        aug_bpc_data.append((t, s, None))
        prev_t = t
        prev_s = s
    t = stop_time_float
    s = 0.0
    for inter in msmsts:
        prop = (inter - prev_t) / (t - prev_t)
        inter_s = prop * s + (1 - prop) * prev_s
        aug_bpc_data.append((inter, inter_s, inter_s))

    if ctype == "txt":
        text = ""
        for (t, s, msms_s) in aug_bpc_data:
            text += "%.3f\t%.0f\n" % (t, s)
    if ctype == "csv":
        text = ""
        for (t, s, msms_s) in aug_bpc_data:
            if msms_s:
                text += "%.3f,%.0f,%.0f\n" % (t, s, msms_s)
            else:
                text += "%.3f,%.0f,\n" % (t, s)
    if ctype == "json":
        text = "["
        for (t, s, msms_s) in aug_bpc_data[:-1]:
            if msms_s:
                text += "[%.3f,%.0f,1],\n" % (t, s)
            else:
                text += "[%.3f,%.0f,0],\n" % (t, s)
        (t, s, msms_s) = aug_bpc_data[-1]
        if msms_s:
            text += "[%.3f,%.0f,1]]\n" % (t, s)
        else:
            text += "[%.3f,%.0f,0]]\n" % (t, s)
    if ctype not in ["txt", "csv", "json"]:
        text = bpc_template(dataset, aug_bpc_data, ctype, filename)
    return text


def do_scan_number_match(filename, scan, sequence, ctype):
    scan = int(scan)
    file = confirm_and_extract(filename)
    the_time = file.scan_list()[scan - 1][0]
    return do_match(filename, the_time, sequence, ctype)


def do_match(dataset, filename, scan_time, sequence, ctype):
    file = confirm_and_extract(dataset, filename)
    the_filters = file.filters()
    (first, last) = file.scan_range()
    scan_num = file.scanForTime(float(scan_time))
    if the_filters[scan_num - 1][1].find("FTMS") > -1:
        mzScan = [(mz, s) for (mz, s, z, n) in file.lscan(scan_num)]
    else:
        mzScan = file.scan(float(scan_time))

    if scan_num == last:
        next_time = scan_time
        prev_time = file.timeForScan(scan_num - 1)
    else:
        next_time = file.timeForScan(scan_num + 1)
        if scan_num == first:
            prev_time = scan_time
        else:
            prev_time = file.timeForScan(scan_num - 1)
    ric_mass = 0.0
    charge = 0
    if the_filters[scan_num - 1][1].find("Full ms ") > -1:
        print "I think this is an MS1..."
        ric_mass = 0.0
    else:
        print "I think this is an MS2..."
        try:
            (ric_mass, charge) = file.scanPrecursor(scan_num)
        except:
            ric_mass = 0.0
            charge = 0
    print "RIC MASS:", ric_mass
    next_ms1 = -1.0
    ms1_nscan = scan_num
    while ms1_nscan != last:
        ms1_nscan += 1
        if the_filters[ms1_nscan - 1][1].find("Full ms ") > -1:
            next_ms1 = file.timeForScan(ms1_nscan)
            break
    last_ms1 = -1.0
    ms1_lscan = scan_num
    while ms1_lscan != first:
        ms1_lscan -= 1
        if the_filters[ms1_lscan - 1][1].find("Full ms ") > -1:
            last_ms1 = file.timeForScan(ms1_lscan)
            break
    scan_type = "cid"
    if the_filters[scan_num - 1][1].find("@etd") > -1:
        scan_type = "etd"
    if the_filters[scan_num - 1][1].find("@hcd") > -1:
        scan_type = "hcd"
    if the_filters[scan_num - 1][1].find("@cid") > -1:
        scan_type = "cid"
    text = match_template(dataset, scan_time, mzScan, sequence, ctype, filename, last_ms1, prev_time, ric_mass,
                          next_time, next_ms1, scan_type, charge, scan_num, ric_mass)
    return text


def do_overlay(dataset, filename, scan_time, sequence, ctype, format):
    file = confirm_and_extract(dataset, filename)
    the_filters = file.filters()
    (first, last) = file.scan_range()
    scan_num = file.scanForTime(float(scan_time))
    if the_filters[scan_num - 1][1].find("FTMS") > -1:
        mzScan = [(mz, s) for (mz, s, z, n) in file.lscan(scan_num)]
    else:
        mzScan = file.scan(float(scan_time))

    if scan_num == last:
        next_time = scan_time
        prev_time = file.timeForScan(scan_num - 1)
    else:
        next_time = file.timeForScan(scan_num + 1)
        if scan_num == first:
            prev_time = scan_time
        else:
            prev_time = file.timeForScan(scan_num - 1)
    charge = 0
    if the_filters[scan_num - 1][1].find("Full ms ") > -1:
        ric_mass = 0.0
    else:
        try:
            (ric_mass, charge) = file.scanPrecursor(scan_num)
        except:
            ric_mass = 0.0
            charge = 0
    next_ms1 = -1.0
    ms1_nscan = scan_num
    while ms1_nscan != last:
        ms1_nscan += 1
        if the_filters[ms1_nscan - 1][1].find("Full ms ") > -1:
            next_ms1 = file.timeForScan(ms1_nscan)
            break
    last_ms1 = -1.0
    ms1_lscan = scan_num
    while ms1_lscan != first:
        ms1_lscan -= 1
        if the_filters[ms1_lscan - 1][1].find("Full ms ") > -1:
            last_ms1 = file.timeForScan(ms1_lscan)
            break
    scan_type = "cid"
    if the_filters[scan_num - 1][1].find("@etd") > -1:
        scan_type = "etd"
    if the_filters[scan_num - 1][1].find("@hcd") > -1:
        scan_type = "hcd"
    if the_filters[scan_num - 1][1].find("@cid") > -1:
        scan_type = "cid"
    text = overlay_template(dataset, scan_time, mzScan, sequence, ctype, filename, last_ms1, prev_time, ric_mass,
                            next_time, next_ms1, scan_type, charge, scan_num, ric_mass, format)
    return text


# Trailing slash in necessary for relative path of img, css and js to work...
@route('/<dataset>/')
def files_request(dataset):
    response.content_type = mimetypes["html"]
    files_template = Template(filename='views/files.html')
    listing = ""
    for f in files(dataset):
        listing += '<a href="files/%s">%s</a><br>' % (f, f)
    text = files_template.render(listing=listing)
    return text


@route('/<dataset>/files/<filename>/scans')
def scanlist_request(dataset, filename):
    response.content_type = mimetypes["html"]
    scans_template = Template(filename='views/scans.html')
    file = confirm_and_extract(dataset, filename)
    scans = file.scan_list(0.0, 1000000.0, 0.0, 1000000.0)
    listing = ""
    for (t, mz) in scans:
        shft = shf(t)
        if mz > 0.0:
            listing += '<a href="%s">%s</a>(<a href="%s">%s</a>)<br>' % (scan_request_url(dataset, filename, t), shft,
                                                                         ric_request_url(dataset, filename, t - 3,
                                                                                         t + 3, mz - 0.05, mz + 0.05),
                                                                         shf(mz))
        else:
            listing += '<a href="%s">%s</a><br>' % (scan_request_url(dataset, filename, t), shft)
    text = scans_template.render(listing=listing)
    return text


@route('/<dataset>/files/<filename>/scans/<scan_time>/header')
def header_request(dataset, filename, scan_time):
    mtype = "txt"  # To simplify url pattern... Also, it makes sense
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_header_request(dataset, filename, scan_time, mtype)


@route('/<dataset>/files/<filename>/scans/<scan_time>/filter')
def filter_request(dataset, filename, scan_time):
    mtype = "txt"  # To simplify url pattern... Also, it makes sense
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_filter_request(dataset, filename, scan_time, mtype)


@route('/<dataset>/files/<filename>/scans/<scan_time>/match/<sequence>')
def match(dataset, filename, scan_time, sequence):
    (sequence, mtype) = untype(sequence)
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_match(dataset, filename, scan_time, sequence, mtype)


@route('/<dataset>/files/<filename>/scans/<scan_time>/overlay/<sequence>')
def full_overlay(dataset, filename, scan_time, sequence):
    (sequence, mtype) = untype(sequence)
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_overlay(dataset, filename, scan_time, sequence, mtype, format="full")


@route('/<dataset>/files/<filename>/scans/scan<scan_number>')
def scan_number_request(dataset, filename, scan_number):
    (scan_number, mtype) = untype(scan_number)
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_scan_number_request(dataset, filename, scan_number, mtype)


@route('/<dataset>/files/<filename>/scans/<scan_time>')
def scan_request(dataset, filename, scan_time):
    (scan_time, mtype) = untype(scan_time)
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_scan_request(dataset, filename, scan_time, mtype)


@route('/<dataset>/files/<filename>/mic/<start_time>-<stop_time>/<start_mz>-<stop_mz>')
def mic_request(dataset, filename, start_time, stop_time, start_mz, stop_mz):
    (stop_mz, mtype) = untype(stop_mz)
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_mic_request(dataset, filename, start_time, stop_time, start_mz, stop_mz, mtype)


@route('/<dataset>/pic/<mz>/<rt>')
def jpic_request(dataset, mz, rt):
    response.content_type = mimetypes["json"]
    return do_jpic_request(dataset, mz, rt)


@route('/<dataset>/xic/<mz>/<rt>')
def pic_request(dataset, mz, rt):
    response.content_type = mimetypes["html"]
    return do_pic_request(dataset, mz, rt)


@route('/<dataset>/files/<filename>/xic/<start_time>-<stop_time>/<start_mz>-<stop_mz>')
def xic_request(dataset, filename, start_time, stop_time, start_mz, stop_mz):
    (stop_mz, mtype) = untype(stop_mz)
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_mic_request(dataset, filename, start_time, stop_time, start_mz, stop_mz, mtype)


@route('/<dataset>/files/<filename>')
def bpc_request(dataset, filename):
    print "filename is", filename
    (filename, mtype) = untype(filename)
    print "now is: ", filename, mtype
    if mtype in mimetypes:
        response.content_type = mimetypes[mtype]
    return do_bpc_request(dataset, filename, mtype)


# Static Routes
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


if __name__ == '__main__':
    run(host='0.0.0.0', port=8888, debug=True)
