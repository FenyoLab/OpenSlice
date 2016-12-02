
function initialize(tolerance_by_source){

    function url_domain(data) {
        var a = document.createElement('a');
        a.href = data;
        return a.hostname;
    }

    var ms1target = null;
    var ms1time = null ;

    var animation = new JS3D("animation") ;

    function l10(val) {
        return Math.log(val) / Math.LN10;
    }

    function e10(val){
        return Math.pow(10,val) ;
    }


    var min_mz = 300.0 ;
    var max_mz = 2000.0 ;

    var lmin = l10(min_mz) ;
    var lmax = l10(max_mz) ;

    var ppm = 1.0e-6 ;

    var hw = l10(1.0 + 20.0*ppm) ;

    var N = ((lmax-lmin)/hw).toFixed() ;

    function smz(mz){
        var lmz = l10(mz)
        if(mz < min_mz || mz > max_mz){
            return 0 ;
        }else{
            var y = ((lmz - lmin) / hw).toFixed() ;
            if(y == 0){
                return 1 ;
            }
            if(y == N){
                return N-1 ;
            }
            return y ;
        }
    }

    a2mod = { "c" : 57.021464, "m": 15.994915, "s": 79.966331 , "t": 79.966331, "y": 79.966331, "p": 15.9949 } ;

    function prot_change(){
        var candidate = d3.select("#protein").property("value") ;
        d3.select("#proteotypic").selectAll("option").remove() ;
        d3.text("http://research.ionomix.com/ptpatlas/"+candidate,function(error,pretext){
            if(error){
                return ;
            }
            ptprows = d3.tsv.parseRows(pretext) ;
            d3.select("#proteotypic").selectAll("option").data(ptprows).enter()
                .append("option").text(function(d) {if(d.length > 1){ return d[1] + " " + d[0];}else{return d[0]} }).property("value",function(d){return d[0]});
        }) ;
    }

    function prot2peps(e){
        e.which = e.which || e.keyCode;
        if(e.which == 13) {
            prot_change() ;
            return false ;
        } else {
            return true ;
        }

    }

    function pset2prots(){
        var pset = d3.select("#proteinset").property("value") ;
        d3.select("#proteinset_details").selectAll("option").remove() ;
        d3.text("../proteinset/"+pset,function(error,pretext){
            if(error){
                return ;
            }
            psetrows = d3.tsv.parseRows(pretext) ;
            d3.select("#proteinset_details").selectAll("option").data(psetrows).enter()
                .append("option").text(function(d) {if(d.length > 1){return d[0] + " " + d[1];}else{return d[0]} }).property("value",function(d){return d[0]});
        }) ;
    }

    function protentry2prot(){
        d3.select("#protein").property("value",d3.select("#proteinset_details").property("value")) ;
        prot_change() ;
    }

    function proteo2peps(){
        d3.select("#peptide").property("value",d3.select("#proteotypic").property("value")) ;
        pep2mass() ;
        load_em_all() ;
    }

    function nist_pep2mass(input_seq){
        var out_varmods = [] ;
        var varmods = [] ;
        var aa_pos = 0 ;
        var clean_seq = "" ;
        for(var i=0; i<input_seq.length; i++){
            aa_pos += 1 ;
            var mod_mass = 0.0 ;
            if(input_seq.charAt(i)=="["){
                var ii=i+1;
                while(ii < input_seq.length){
                    if(input_seq.charAt(ii)=="]"){
                        mod_mass = parseFloat(input_seq.substring(i+1,ii)) ;
                        i = ii + 1;
                        break ;
                    } else {
                        ii += 1 ;
                    }
                }
            }
            var c = input_seq.charAt(i) ;
            clean_seq = clean_seq + c.toUpperCase() ;
            if(c == c.toLowerCase()){
                varmods.push(new VariableModification(aa_pos,a2mod[c]+mod_mass,AminoAcid.get(c.toUpperCase())));
                out_varmods.push({"index": aa_pos,"modMass": a2mod[c]+mod_mass, "aminoAcid": c.toUpperCase()}) ;
            }
            if( (c != c.toLowerCase()) && (mod_mass != 0.0) ){
                varmods.push(new VariableModification(aa_pos,mod_mass,AminoAcid.get(c.toUpperCase())));
                out_varmods.push({"index": aa_pos,"modMass": mod_mass, "aminoAcid": c.toUpperCase()}) ;
            }
        }

        var ignored = new Peptide( clean_seq , null, varmods, 0.0, 0.0 ) ;
        var int_z = parseInt(d3.select("#spinner").property("value")) ;
        var z = parseFloat(d3.select("#spinner").property("value")) ;
        // Note that Ion.MASS_H is not the same as Ion.MASS_PROTON!!! (Thanks to Robert Chalkley for catching this one)...
        var new_mz =  (Peptide.getNeutralMassMono()+ z * Ion.MASS_PROTON) / z  ;
        return {"mz": new_mz, "varmods": out_varmods, "seq": clean_seq, "z": int_z} ;
    }



    function data_pep2mass(){
        var out_varmods = [] ;
        var varmods = [] ;
        var input_seq = d3.select("#peptide").property("value") ;
        var aa_pos = 0 ;
        var clean_seq = "" ;
        for(var i=0; i<input_seq.length; i++){
            aa_pos += 1 ;
            var mod_mass = 0.0 ;
            if(input_seq.charAt(i)=="["){
                var ii=i+1;
                while(ii < input_seq.length){
                    if(input_seq.charAt(ii)=="]"){
                        mod_mass = parseFloat(input_seq.substring(i+1,ii)) ;
                        i = ii + 1;
                        break ;
                    } else {
                        ii += 1 ;
                    }
                }
            }
            var c = input_seq.charAt(i) ;
            clean_seq = clean_seq + c.toUpperCase() ;
            if(c == c.toLowerCase()){
                varmods.push(new VariableModification(aa_pos,a2mod[c]+mod_mass,AminoAcid.get(c.toUpperCase())));
                out_varmods.push({"index": aa_pos,"modMass": a2mod[c]+mod_mass, "aminoAcid": c.toUpperCase()}) ;
            }
            if( (c != c.toLowerCase()) && (mod_mass != 0.0) ){
                varmods.push(new VariableModification(aa_pos,mod_mass,AminoAcid.get(c.toUpperCase())));
                out_varmods.push({"index": aa_pos,"modMass": mod_mass, "aminoAcid": c.toUpperCase()}) ;
            }
        }

        var ignored = new Peptide( clean_seq , null, varmods, 0.0, 0.0 ) ;
        var int_z = parseInt(d3.select("#spinner").property("value")) ;
        var z = parseFloat(d3.select("#spinner").property("value")) ;
        // Note that Ion.MASS_H is not the same as Ion.MASS_PROTON!!! (Thanks to Robert Chalkley for catching this one)...
        var new_mz =  (Peptide.getNeutralMassMono()+ z * Ion.MASS_PROTON) / z  ;
        return {"mz": new_mz, "varmods": out_varmods, "seq": clean_seq, "z": int_z} ;
    }

    function pep2mass_non_event(){
        d3.select("#mz").property("value",data_pep2mass().mz.toFixed(5)) ;
        load_em_all() ;
    }

    function pep2mass(e){
        d3.select("#mz").property("value",data_pep2mass().mz.toFixed(5)) ;
        if(e){
            e.which = e.which || e.keyCode;
            if(e.which == 13) {
                load_em_all() ;
                return false ;
            } else {
                return true ;
            }
        }
    }

    function return_launch(e){
        if(e){
            e.which = e.which || e.keyCode;
            if(e.which == 13) {
                load_em_all() ;
                return false ;
            } else {
                return true ;
            }
        }
    }

    function score(targetURL,input_seq,tolerance,callback){
        d3.text(targetURL,"text/csv",function(the_text){
            var peaks = [] ;
            var total = 0 ;
            var score = 0.0 ;

            rows = d3.csv.parseRows(the_text) ;
            rows.shift() ; // precursor,charge,ignore <-- BTW, should be precursor,ignore,charge (to match the rest of the file)
            rows.forEach(function(r){
                var intensity = +r[1] ;
                var charge = +r[2] ;
                if(tolerance < 0.1){
                    if(charge > 0){
                        total += intensity ;
                    }
                }else{
                    total += intensity ;
                }
                peaks.push({"mz": +r[0], "intensity": intensity, "charge": charge, "dirty": false}) ;
            }) ;

            var out_varmods = [] ;
            var varmods = [] ;
            var aa_pos = 0 ;
            var clean_seq = "" ;
            for(var i=0; i<input_seq.length; i++){
                aa_pos += 1 ;
                var mod_mass = 0.0 ;
                if(input_seq.charAt(i)=="["){
                    var ii=i+1;
                    while(ii < input_seq.length){
                        if(input_seq.charAt(ii)=="]"){
                            mod_mass = parseFloat(input_seq.substring(i+1,ii)) ;
                            i = ii + 1;
                            break ;
                        } else {
                            ii += 1 ;
                        }
                    }
                }
                var c = input_seq.charAt(i) ;
                clean_seq = clean_seq + c.toUpperCase() ;
                if(c == c.toLowerCase()){
                    varmods.push(new VariableModification(aa_pos,a2mod[c]+mod_mass,AminoAcid.get(c.toUpperCase())));
                    out_varmods.push({"index": aa_pos,"modMass": a2mod[c]+mod_mass, "aminoAcid": c.toUpperCase()}) ;
                }
                if( (c != c.toLowerCase()) && (mod_mass != 0.0) ){
                    varmods.push(new VariableModification(aa_pos,mod_mass,AminoAcid.get(c.toUpperCase())));
                    out_varmods.push({"index": aa_pos,"modMass": mod_mass, "aminoAcid": c.toUpperCase()}) ;
                }
            }

            var ignored = new Peptide( clean_seq , null, varmods, 0.0, 0.0 ) ;
            var int_z = parseInt("2") ;
            var z = parseFloat("2") ;
            // Note that Ion.MASS_H is not the same as Ion.MASS_PROTON!!! (Thanks to Robert Chalkley for catching this one)...
            var new_mz =  (Peptide.getNeutralMassMono()+ z * Ion.MASS_PROTON) / z  ;

            var ions = [] ;
            for(i = 1 ; i < input_seq.length ; i++){
                for(z = 1 ; z <= 3 ; z++){
                    ions.push({"mz": Ion.getSeriesIon({"charge": z, "type": "y"}, input_seq, i , "mono").mz, "charge": z}) ;
                    ions.push({"mz": Ion.getSeriesIon({"charge": z, "type": "b"}, input_seq, i , "mono").mz, "charge": z}) ;
                }
            }

            ions.forEach(function(tion){
                for(var i = 0 ; i < peaks.length ; i++){
                    var p = peaks[i] ;
                    if(p.dirty){
                        continue ;
                    }
                    if(Math.abs(tion.mz-p.mz)<=tolerance){
                        if(tolerance < 0.1){
                            if(p.charge != tion.charge){
                                continue ;
                            }
                        }
                        p.dirty = true ;
                    }
                }
            })

            var explained = 0 ;
            for(i = 0 ; i < peaks.length ; i++){
                p = peaks[i] ;
                if(p.dirty){
                    explained += p.intensity ;
                }
            }
            score = explained / total ;
            callback(score) ;
        });
    }

    var file_times = d3.map() ;
    var xics = d3.map() ;


    function pic(target_mz,callback){

        var thread_times = function(){
            d3.json("pic/" + smz(target_mz),function(error,raw_xics){
                if(error){
                    callback() ;
                } else {
                    all_files.forEach(function(fn){
                        var the_vals = raw_xics[fn] ;
                        var ts = file_times.get(fn) ;
                        answer = [] ;
                        time_offset = 0 ;
                        for(var i = 0; i < the_vals.length ; i++){
                            var cast = +the_vals[i] ;
                            if(cast < 0){
                                answer.push({"time": ts[time_offset], "intensity": 0}) ;
                                var last_time_offset = time_offset ;
                                time_offset = -1*(cast+1) ;
                                if(time_offset-1 > last_time_offset){
                                    answer.push({"time": ts[time_offset-1], "intensity": 0}) ;
                                }
                            } else {
                                answer.push({"time": ts[time_offset], "intensity": cast}) ;
                                time_offset += 1 ;
                            }
                        }
                        xics.set(fn,{"file": fn, "mz": target_mz, "values": answer})
                    }) ;
                    callback() ;
                }
            }) ;
        } ;

        thread_times() ;
    }


    function xic(target_file,target_mz,callback){

        var thread_times = function(){
            d3.text("files/" + target_file + "/xics/" + smz(target_mz),function(error,the_text){
                if(error){
                    xics.set(target_file,{"file": target_file, "mz": target_mz, "values": []})
                    callback() ;
                } else {
                    var the_vals = d3.merge(d3.csv.parseRows(the_text)) ;
                    var ts = file_times.get(target_file) ;
                    answer = [] ;
                    time_offset = 0 ;
                    for(var i = 0; i < the_vals.length ; i++){
                        var cast = +the_vals[i] ;
                        if(cast < 0){
                            answer.push({"time": ts[time_offset], "intensity": 0}) ;
                            var last_time_offset = time_offset ;
                            time_offset = -1*(cast+1) ;
                            if(time_offset-1 > last_time_offset){
                                answer.push({"time": ts[time_offset-1], "intensity": 0}) ;
                            }
                        } else {
                            answer.push({"time": ts[time_offset], "intensity": cast}) ;
                            time_offset += 1 ;
                        }
                    }
                    xics.set(target_file,{"file": target_file, "mz": target_mz, "values": answer})
                    callback() ;
                }
            }) ;
        } ;

        thread_times() ;
    }

    var margin = {top: 20, right: 40, bottom: 100, left: 70},
        margin2 = {top: 450, right: 40, bottom: 30, left: 70},
        width = 940 - margin.left - margin.right,
        height = 500 - margin.top - margin.bottom,
        height2 = 500 - margin2.top - margin2.bottom;

    var x = d3.scale.linear().range([0, width]),
        x2 = d3.scale.linear().range([0, width]),
        y = d3.scale.linear().range([height, 0]),
        y2 = d3.scale.linear().range([height2, 0]) ;


    color = function(fname){
        return published.colors[published[fname][active_category]] ;
    }

    var xAxis  = d3.svg.axis().scale(x).orient("bottom"),
        xAxis2 = d3.svg.axis().scale(x2).orient("bottom") ;
    var yAxis = d3.svg.axis().scale(y).orient("left");

    var line = d3.svg.line()
        .interpolate("linear")
        .x(function(d) { return x(d.time); })
        .y(function(d) { return y(d.intensity); });

    var line2 = d3.svg.line()
        .interpolate("linear")
        .x(function(d) { return x2(d.time); })
        .y(function(d) { return y2(d.intensity); });

    var ms1_svg = d3.select("#ms1_holder").append("svg") ;

    var svg = d3.select("body").append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom) ;

    svg.append("defs").append("clipPath")
        .attr("id", "clip")
        .append("rect")
        .attr("width", width)
        .attr("height", height);

    var focus = svg.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")").attr("id","ze_focus");
    var context = svg.append("g")
        .attr("transform", "translate(" + margin2.left + "," + margin2.top + ")");

    focus.append("g")
        .attr("id","xaxis")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis)
        .append("text")
        .attr("x",width-margin.right-20)
        .attr("y", 30)
        .text("Retention Time (mins)");

    focus.append("g")
        .attr("id","yaxis")
        .attr("class", "y axis")
        .call(yAxis)
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("Intensity");

    context.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height2 + ")")
        .call(xAxis2);

    function brushed() {
        x.domain(brush.empty() ? x2.domain() : brush.extent());
        display() ;
    }

    var brush = undefined ;

    var target = svg.append("g")
        .attr("class", "target")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    target.append("line").attr("display", "none")
        .attr("class", "targetx")
        .attr("y1", y(0) - 6)
        .attr("y2", y(0) + 6);

    target.append("circle").attr("display", "none")
        .attr("class", "targetcircle").attr("r", 3.5) ;

    target.append("text").attr("id","targetlabel").text("").attr("display", "none");

    target.on("click",function(d){
        var altshift = false ;
        if (d3.event.shiftKey) {
            altshift = true ;
                if(ms1target){
                    var survives = published[ms1target.file][active_category] ;
                    categories[active_category].forEach(function(val){
                        if(val != survives){
                            multifilter[val] = -1 ;
                        }
                    }) ;
                }
        }
        if (d3.event.altKey) {
            altshift = true ;
                if(ms1target){
                    var survives = published[ms1target.file][active_category] ;
                    published.variables.forEach(function(v){
                        categories[v.name].forEach(function(val){
                            if(published[ms1target.file][v.name] != val){
                                multifilter[val] = -1 ;
                            }
                        })
                    }) ;                 
                }
        }
        if(altshift){
            transdisplay()
        }else{
            target_domain = url_domain(document.URL) ;
            dataset = document.URL.split("/").slice(-2,-1)[0] ;
            https = document.URL.slice(0, document.URL.indexOf(":")) ;
            window.open(https + "://" + target_domain + ":8888/" + dataset + "/files/" + ms1target.file +  "/scans/" + ms1time) ;
        }
    }) ;

    svg.append("text").attr("id","ms2label").text("").attr("display", "none");

    target.append("rect")
        .attr("class", "targetoverlay")
        .attr("width", width)
        .attr("height", height)
        .on("mousemove", mousemove)
        .on("mouseleave", hidetarget);

    function hidetarget(){
        target.select(".targetcircle").attr("display", "none");
        target.select(".targetx").attr("display", "none");
        target.select("#targetlabel").attr("display","none")
    }

    function mousemove() {

        var in_x = d3.mouse(this)[0] ;
        var in_y = d3.mouse(this)[1] ;

        var the_x = x.invert(d3.mouse(this)[0]) ;
        var the_y = y.invert(d3.mouse(this)[1]) ;

        var bisect = d3.bisector(function(d) { return d.time; }).right;
        var candidates = [] ;
        xics.forEach(function(fname,xc){
            if(multifiltered(fname)){
                return ;
            }
            var offset = bisect(xc.values,the_x) ;
            if(offset < xc.values.length){
                var dp = xc.values[offset] ;
                candidates.push({"time": dp.time, "intensity": dp.intensity, "file": fname}) ;
            }
            if(offset > 0){
                var dp = xc.values[offset-1] ;
                candidates.push({"time": dp.time, "intensity": dp.intensity, "file": fname}) ;
            }
        });

        if(candidates.length > 0){
            var closest = candidates[0] ;
            var cc_x = x(closest.time) ;
            var cc_y = y(closest.intensity) ;
            var closest_d_c = (in_x - cc_x)*(in_x - cc_x) + (in_y - cc_y)*(in_y - cc_y)
            candidates.forEach(function(cand){
                var c_x = x(cand.time) ;
                var c_y = y(cand.intensity) ;
                var d_c = (in_x - c_x)*(in_x - c_x) + (in_y - c_y)*(in_y - c_y) ;
                if(d_c < closest_d_c){
                    closest_d_c = d_c ;
                    closest = cand ;
                }
            }) ;
            c_x = x(closest.time) ;
            c_y = y(closest.intensity) ;
            c_name = closest.file ;
            ms1target = closest ;
            ms1time = closest.time ;

            target.select(".targetcircle").attr("transform", "translate(" + c_x + "," + c_y + ")");
            target.select(".targetx").attr("transform", "translate(" + c_x + ",0)");
            target.select(".targetcircle").attr("display", "inline") ;
            target.select("#targetlabel").text(published[closest.file].label).attr("fill","black").attr("transform", "translate(" + (c_x + 10) + "," + (c_y +5) + ")").attr("display", "inline") ;
            target.select(".targetx").attr("display", "inline");
        } else {
            ms1target = null ;
            target.select(".targetcircle").attr("display", "none");
            target.select("#targetlabel").attr("display","none") ;
            target.select(".targetx").attr("display", "none");
        }
    }

    var score_threshold = 0.2 ;

    function score_all_ms2s(){
        if(d3.select("#offscore").property("checked") == true){
            score_timeout_var = setTimeout(score_all_ms2s,50) ;
        } else {
            var the_filter = null ;
            if(d3.select("#allscore").property("checked") == true){
                the_filter = function(d){return d[9] == -1.0 ;} ;
            } else {
                the_filter = function(d){return d[1] > x.domain()[0] && d[1] < x.domain()[1] && !d[7] && d[9] == -1.0 ;} ;
            }

            var candidates = [] ;
            if(ms2pos){
                candidates = ms2pos.filter(the_filter) ;
            }

            if(candidates.length){
                d = candidates[0] ;
                var target_url = "files/"+d[0]+"/scans/"+d[5] ;
                var target_seq = hidden_seq ;
                score(target_url,target_seq,tolerance_by_source(d[0]),function(the_score){
                    d[9] = the_score ;
                    if(the_score >= score_threshold){
                        d3.select("#file_" + d[0] +"_"+d[1].replace(".","_")).attr("height",8) ;
                    }else{
                        d3.select("#file_" + d[0] +"_"+d[1].replace(".","_")).attr("height",4) ;
                    }
                    inject_report() ;
                    score_timeout_var = setTimeout(score_all_ms2s,50)
                });
            } else {
                score_timeout_var = setTimeout(score_all_ms2s,50)
            }

        }
        if(ms2pos){
            ms2pos.filter(function(d){return d[7]}).forEach(function(d){
                if(d[7]!==d3.select("#peptide").property("value").toUpperCase()){
                    d3.select("#" + d[0] +"_"+d[1].replace(".","_")).attr("height",4) ;
                }else{
                    d3.select("#" + d[0] +"_"+d[1].replace(".","_")).attr("height",8) ;
                }
            });
        }
    }

    var the_bpc = [] ;
    //This variable is an incredibly wrong hack... Fix this at some point...
    var legit_iframe_load = false ;
    function complete_display(){
        if(!focus.selectAll(".xic").empty()){
            focus.selectAll(".xic").remove() ;
        }
        if(!target.selectAll(".point").empty()){
            target.selectAll(".point").remove() ;
        }
        if(!context.select("#bpc").empty()){
            context.select("#bpc").remove() ;
        }
        if(!context.select("#xbrush").empty()){
            context.select("#xbrush").remove() ;
        }
        var bpc = d3.map() ;
        xics.forEach(function(fname,xc){
            xc.values.forEach(function(v){
                var deci = +(v.time.toFixed(1)) ;
                if(!bpc.has(deci)){
                    bpc.set(deci, v.intensity) ;
                    var before = +((deci-0.1).toFixed(1)) ;
                    var after = +((deci+0.1).toFixed(1)) ;
                    if(!bpc.has(before)){
                        bpc.set(before,0) ;
                    }
                    if(!bpc.has(after)){
                        bpc.set(after,0) ;
                    }
                }else{
                    bpc.set(deci, d3.max([bpc.get(deci), v.intensity])) ;
                    var before = +((deci-0.1).toFixed(1)) ;
                    var after = +((deci+0.1).toFixed(1)) ;
                    if(!bpc.has(before)){
                        bpc.set(before,0) ;
                    }
                    if(!bpc.has(after)){
                        bpc.set(after,0) ;
                    }
                }
            }) ;
        }) ;
        the_bpc = [] ;
        the_bpc_keys = bpc.keys() ;
        the_bpc_keys.sort(function(a,b){return a-b}) ;
        bpc_times = the_bpc_keys.forEach(function(a_key){
            the_bpc.push({"time": a_key, "intensity": bpc.get(a_key)}) ;
        }) ;

        var x_min = d3.min(xics.values(), function(an_xic){ return d3.min(an_xic.values, function(v){ return v.time; });}) ;
        var x_max = d3.max(xics.values(), function(an_xic){ return d3.max(an_xic.values, function(v){ return v.time; });}) ;
        x2.domain([x_min,x_max]) ;

        var y_min = 0 ;
        var y_max = 0
        var y_max_time = 0 ;
        xics.values().forEach(function(an_xic){
            an_xic.values.forEach(function(v){
                if(v.intensity > y_max){
                    y_max = v.intensity ;
                    y_max_time = v.time ;
                }
            })
        })
        y2.domain([y_min , y_max]) ;

        var brush_start = y_max_time - 8 ;
        var brush_stop = y_max_time + 8 ;

        x.domain([brush_start,brush_stop]) ;
        y.domain([0,y_max]) ;

        var area2 = d3.svg.area().interpolate("monotone").x(function(d){return x2(d.time)}).y0(height2).y1(function(d){return y2(d.intensity)}) ;
        context.append("path").attr("id","bpc").attr("d",area2(the_bpc)) ;

        var xic_elems = focus.selectAll(".xic").data(xics.values()) ;
        xic_elems.enter().append("g").attr("class","xic").append("path").attr("class","line").attr("clip-path", "url(#clip)") ;

        function score_blurb(a_score){
            if(a_score < 0){
                return "" ;
            }else{
                return " (" + a_score.toFixed(2) + ")" ;
            }
        }
        var points = target.selectAll(".point")
            .data(ms2pos)
            .enter().append("svg:rect")
            .attr("class","point")
            .attr("id",function(d){return "file_" + d[0]+"_"+d[1].replace(".","_")})
            .attr("fill", function(d, i) { return color(d[0]) })
            .attr("stroke", function(d, i) { if(d[7]){return "gray"}else{return color(d[0])} })
            .attr("stroke-width",function(d,i){if(d[7]){return 2}else{return 1}})
            .attr("x", function(d, i) { return x(d[1]) })
            .attr("y", function(d, i) { return "-10px" })
            .attr("width", function(d, i) { return 6 })
            .attr("height", function(d, i) { return 6 })
            .on("click",function(d){
                d3.event.stopPropagation();


                var info = data_pep2mass() ;
                var nist_info = nist_pep2mass(d[7])

                var altshift = false ;
                if (d3.event.shiftKey) {
                    altshift = true ;
                        var survives = published[d[0]][active_category] ;
                        categories[active_category].forEach(function(val){
                            if(val != survives){
                                multifilter[val] = -1 ;
                            }
                        }) ;
                }
                if (d3.event.altKey) {
                        var survives = published[d[0]][active_category] ;
                        published.variables.forEach(function(v){
                            categories[v.name].forEach(function(val){
                                if(published[ms1target.file][v.name] != val){
                                    multifilter[val] = -1 ;
                                }
                            })
                        }) ;                 
                }
                if(altshift){
                    transdisplay() ;
                    return ;
                }

                var targetscan = "files/"+d[0]+"/scans/"+d[5] ;
                localStorage["scanURL"] = targetscan ;
                localStorage["sequence"] = info.seq ;
                if(d[7]){
                    localStorage["nist_peptide"] = d[7] ;
                    localStorage["nist_URL"] = "../nist/"+d[6] ;
                    localStorage["nist_variableMods"] = "[]" ;
                }else{
                    localStorage["nist_peptide"] = "" ;
                    localStorage["nist_variableMods"] = [] ;
                    localStorage["nist_URL"] = "" ;
                }
                localStorage["scanTime"] = d[1] ;
                localStorage["fileName"] = d[0] ;
                localStorage["charge"] = info.z ;
                localStorage["variableMods"] = JSON.stringify(info.varmods) ;
                localStorage["massError"] = tolerance_by_source(d[0]) ;
                legit_iframe_load = true ;
                window.open("ms2_viewer.html","firstmsmsview","menubar=0,location=0,resizable=0,scrollbars=0,status=0") ;
                window.open("nist_viewer.html","nistview","menubar=0,location=0,resizable=0,scrollbars=0,status=0") ;
                if(d[7]){
                d3.select("#nist_button").style("display","inline") ;
                }else{
                d3.select("#nist_button").style("display","none") ;
                }
            })
            .on("mouseover", function(d){
                if(d[7]){
                    d3.select("#ms2label").text(published[d[0]].label+" [NIST: "+ d[7] + "]" +score_blurb(d[9])).attr("display","inline").attr("fill","black").attr("transform", "translate(" + (x(d[1])+40) + "," + (8) + ")") ; d3.select(this).style("stroke", "black" );
                }else{
                    d3.select("#ms2label").text(published[d[0]].label+score_blurb(d[9])).attr("display","inline").attr("fill","black").attr("transform", "translate(" + (x(d[1])+40) + "," + (8) + ")") ; d3.select(this).style("stroke", "black" );
                }
            })
            .on("mouseout", function(d){
                if(d[7]){
                    d3.select("#ms2label").attr("display","none"); d3.select(this).style("stroke", "gray");
                }else{
                    d3.select("#ms2label").attr("display","none"); d3.select(this).style("stroke", color(d[0]));
                }
            });


        brush = d3.svg.brush()
            .x(x2) ;

        context.append("g")
            .attr("class", "x brush")
            .attr("id","xbrush")
            .call(brush.extent([brush_start,brush_stop]))
            .selectAll("rect")
            .attr("height", height2);

        brush.on("brush", brushed);

        context.select(".x.axis").call(xAxis2) ;
        display() ;
    }

    function display(){

        var grand_max = 0 ;
        var visibles = d3.set() ;
        var ms2_visibles = d3.set() ;

        var grand_maxer = function(v){
            if((v.time > x.domain()[0]) && (v.time < x.domain()[1])){
            }
        } ;

        var zero_filter = function(v){ if( (v.time > x.domain()[0]) && (v.time < x.domain()[1])){return v.intensity;}else{return 0.0}} ;

        var dfilter = function(stuff){
            substuff = [] ;
            stuff.forEach(function(v){
                if( (v.time > x.domain()[0]) && (v.time < x.domain()[1]) ){
                    if(v.intensity > grand_max){
                       grand_max = v.intensity ;
                    }
                    substuff.push(v) ;
                }
            }) ;
            return substuff ;
        } ;

        y.domain([0.0,d3.max(the_bpc,zero_filter)]) ;

        var xic_elems = focus.selectAll(".xic").each(function(d){
            if(multifiltered(d.file)){
                return ;
            }
            visibles.add(d.file) ;
            var subdata = dfilter(d.values) ;
            if(subdata.length == 0){
                return ;
            } else {
                return ;
            }
        }) ;

        y.domain([0.0,grand_max]) ;

        var xic_elems = focus.selectAll(".xic") ;
        xic_elems.select("path").attr("d", function(d){
            if(multifiltered(d.file)){
                return line([]);
            }
            var subdata = dfilter(d.values) ;
            if(subdata.length == 0){
                //Should be the same as above...
                return line([]) ;
            } else {
                return line(subdata);
            }
        }).style("stroke", function(d){return color(d.file);});

        target.selectAll(".point")
            .attr("x", function(d, i) { return x(d[1]) })
            .attr("y", function(d, i) { return "-10px" })
            .attr("display", function(d, i){
                if(multifiltered(d[0])){
                    return "none" ;
                }

                // Count ms2 scans even if they are out of the time-domain...
                ms2_visibles.add(d) ;
                if( (d[1] > x.domain()[0]) && (d[1] < x.domain()[1])){
                    return "inline"
                } else {
                    return "none"
                }
            }
            ) ;

        focus.select(".x.axis").call(xAxis) ;
        focus.select(".y.axis").call(yAxis) ;

        hidetarget() ;

        d3.select("#totals").html("Raw Files: " + visibles.values().length + " MS2 Scans: " + ms2_visibles.values().length) ;
        inject_report() ;
    }

    function multifiltered(fname){
        var retval = false ;
        published.variables.forEach(function(v){
            if(multifilter[published[fname][v.name]] < 0){
                retval = true ;
            }
        }) ;
        return retval ;
    }

    function transdisplay(){

    var grand_max = 0 ;
    var visibles = d3.set() ;
    var ms2_visibles = d3.set() ;

    var grand_maxer = function(v){
        if((v.time > x.domain()[0]) && (v.time < x.domain()[1])){
        }
    };

    var zero_filter = function(v){ if( (v.time > x.domain()[0]) && (v.time < x.domain()[1])){return v.intensity;}else{return 0.0}} ;

    var dfilter = function(stuff){
        substuff = [] ;
        stuff.forEach(function(v){
            if( (v.time > x.domain()[0]) && (v.time < x.domain()[1]) ){
                if(v.intensity > grand_max){
                   grand_max = v.intensity ;
                }
                substuff.push(v) ;
            }
        }) ;
        return substuff ;
    } ;

    y.domain([0.0,d3.max(the_bpc,zero_filter)]) ;

    var xic_elems = focus.selectAll(".xic").each(function(d){
        if(multifiltered(d.file)){
            return line([]);
        }
        visibles.add(d.file) ;
        var subdata = dfilter(d.values) ;
        if(subdata.length == 0){
            return ;
        } else {
            return ;
        }
    }) ;

    y.domain([0.0,grand_max]) ;

    var xic_elems = focus.selectAll(".xic") ;
    xic_elems.select("path").transition().duration(1000).attr("d", function(d){
        if(multifiltered(d.file)){
            return line([]);
        }

        var subdata = dfilter(d.values) ;
        if(subdata.length == 0){
            //Should be the same as above...
            return line([]) ;
        } else {
            return line(subdata);
        }
    }).style("stroke", function(d){return color(d.file);});

    target.selectAll(".point")
        .attr("x", function(d, i) { return x(d[1]) })
        .attr("y", function(d, i) { return "-10px" })
        .attr("fill", function(d, i){ return color(d[0])})
        .attr("stroke", function(d, i) { if(d[7]){return "gray"}else{return color(d[0])} })
        .attr("display", function(d, i){

            if(multifiltered(d[0])){
                return "none" ;
            }

            // Count ms2 scans even if they are out of the time-domain...
            ms2_visibles.add(d) ;
            if( (d[1] > x.domain()[0]) && (d[1] < x.domain()[1])){

                return "inline"
            } else {
                return "none"
            }
        }) ;

    focus.select(".x.axis").call(xAxis) ;
    focus.select(".y.axis").call(yAxis) ;

    hidetarget() ;

    d3.select("#totals").html("Raw Files: " + visibles.values().length + " MS2 Scans: " + ms2_visibles.values().length) ;
    inject_report() ;
}

var multifilter = {} ;

var active_category = null ;

function change_active_category(new_cat){
    if(new_cat != active_category){
        active_category = new_cat ;
        transdisplay() ;
    }
}

function clear_variable(v){
    categories[v].forEach(function(cat){
        multifilter[cat] = 1 ;
    })
    transdisplay() ;
}

function clear_multifilter(){
    published.variables.forEach(function(v){
        categories[v.name].forEach(function(cat){
            multifilter[cat] = 1 ;
        })
    })
    transdisplay() ;
}

function flip_category(value){
    multifilter[value] *= -1 ;
    transdisplay() ;
}

function multifilter_strikethrough(cat){
    if(multifilter[cat] > 0){
        return "" ;
    } else {
        return " text-decoration: line-through; "
    }
}

function inject_report(){
    xic_map = d3.map() ;
    scan_map = d3.map() ;

    if(xics.keys().length < all_files.length){
        return ;
    }

    xics.keys().forEach(function(k){
        var catkey = published[k][active_category] ;
        if(!xic_map.has(catkey)){
            xic_map.set(catkey,0)
        }
        xic_map.set(catkey,xic_map.get(catkey)+1)
    }) ;
    ms2pos.forEach(function(r){
        var catkey = published[r[0]][active_category] ;
        if(!scan_map.has(catkey)){
            scan_map.set(catkey,[0,0,0])
        }
        var seen = scan_map.get(catkey)[0]
        var scored = scan_map.get(catkey)[1]
        var hit = scan_map.get(catkey)[2]
        seen += 1 ;
        if(r[9]>=0){
            scored += 1 ;
            if(r[9]>=score_threshold){
                hit += 1 ;
            }
        }
        scan_map.set(catkey,[seen,scored,hit])
    }) ;
    report_html = "" ;

    published.variables.forEach(function(entry){
        if(entry.name == active_category){
            report_html = report_html + "<tr><td onclick='viewer.clear_variable(\"" + entry.name + "\")' style='cursor: crosshair; text-align: right; font-weight: bold; font-style: sans-serif'>" + entry.name + ":</td><td>"
            entry.values.forEach(function(k){
                   report_html = report_html +
                    "<span onclick='viewer.flip_category(\"" + k + "\")' style='cursor: crosshair; text-align: right;font-style: sans-serif; color: " + published.colors[k] + ";" + multifilter_strikethrough(k) + "'>" + k + "</span>&nbsp;" 
            })
            report_html = report_html + "</td></tr>"
        } else {
            report_html = report_html + "<tr><td onclick='viewer.change_active_category(\"" + entry.name + "\")' style='cursor: crosshair; text-align: right; font-style: sans-serif'>" + entry.name + ":</td><td>"
            entry.values.forEach(function(k){
                   report_html = report_html +
                    "<span onclick='viewer.flip_category(\"" + k + "\")' style='cursor: crosshair; text-align: right;font-style: sans-serif;"+  multifilter_strikethrough(k) +"'>" + k + "</span>&nbsp;" 
            })
            report_html = report_html + "</td></tr>"
        }
    })

    d3.select("#report_table").html(
        report_html
    )
}



    function expose_ms2_frame(){
        if(legit_iframe_load){
            d3.select('#ms2_iframe_holder').style('opacity',0.0) ;
            d3.select('#ms2_iframe_holder').style('display','inline') ;
            d3.select('#ms2_iframe_holder').transition().duration(500).style('opacity',1.0) ;
        }
    }
    function expose_nist_frame(){
        if(legit_iframe_load){
            d3.select('#nist_iframe_holder').style('display','inline') ;
        }
    }
    function renderloop() {
        animation.matrix = animation.rotateX(0.05);
        animation.matrix = animation.rotateZ(-0.05);
        animation.matrix = animation.rotateY(0.05);
        animation.paint();
        if(file_times.keys().length == all_files.length){
            return false ;
        } else {
            return true ;
        }
    }

    function animLoop( render ) {
        var running, lastFrame = +new Date,
            raf = window.requestAnimationFrame || window.mozRequestAnimationFrame    ||
                window.webkitRequestAnimationFrame ||
                window.msRequestAnimationFrame     ||
                window.oRequestAnimationFrame;
        function loop( now ) {
            // stop the loop if render returned false
            if ( running !== false ) {
                raf( loop );
                var deltaT = now - lastFrame;
                if ( deltaT < 100000 ) {
                    running = render(  );
                }
                lastFrame = now;
            }
        }
        loop( lastFrame );
    }

    var pointmap = d3.map()
    function initialize_times(file_list,finisher){
        var launched = 0 ;
        file_list.forEach(function(fname){
            var ax1 = d3.shuffle([0,1])[0] ;
            var ax2 = d3.shuffle([0,1])[0] ;
            var ax3 = Math.random()*1 ;
            axs = [ax1,ax2,ax3] ;
            axs = d3.shuffle(axs) ;
            animation.matrix = identity() ;
            animation.matrix = animation.translate( document.getElementById("animation").offsetWidth/2 -380,200,-400);
            m = animation.scale(300,300,300,identity())
            pointmap.set(fname,animation.addPoint(axs[0],axs[1],-1*axs[2],fname,m)) ;
        }) ;
        animLoop(renderloop) ;

        file_list.forEach(function(fname){
            pointmap.get(fname).div.style.color = color(fname)  ;
        }) ;

        d3.json("xic_times",function(error,json){
            file_times = d3.map(json) ;
            finisher() ;
        }) ;
    }


    var all_files = undefined ;

    var published = undefined ;
    var categories = {} ;

    var ms2pos = undefined ;
    var hidden_mz = undefined ;
    var hidden_seq = undefined ;
    var hidden_charge = undefined ;
    function load_em_all(mz_query){
        if(!mz_query){
            mz_query = parseFloat(d3.select("#mz").property("value")).toFixed(5) ;
        }
        hidden_mz = mz_query ;
        hidden_seq = d3.select("#peptide").property("value") ;
        hidden_charge = d3.select("#spinner").property("value") ;

        d3.select("#launch_image").attr("src","img/dice_rotate.gif") ;
	d3.select("#diceintitle").attr("src","img/dice_icon_small_rotating.gif") ;
        var after_initialize_times = function(){
            d3.select("#waitmessage").html("") ;
	    d3.select("#animation").html("<table height='100%' width='100%'><tr width='100%' height='100%'><td width='100%' height='100%' align='center'><span class='emergency'>&lt;&lt; extracting <strong>G</strong>lobal <strong>I</strong>on <strong>C</strong>hromatogram &gt;&gt;</span> </td></tr></table>") ;
            if(d3.select("#animation").style("display") !== "none"){
              d3.select("#animation").transition().duration(1000).style('opacity',0.5) ;
            }else{
              d3.select("#animation").style("display","inline") ;
              d3.select("#animation").transition().duration(1000).style('opacity',0.5) ;
            }

            xics = d3.map() ;
            pic(mz_query,function(){
                ms2pos = [] ;
                d3.text("precursors/"+smz(mz_query),function(error,pretext){
                    if(error){
                        ms2pos = [] ;
                    } else {
                        ms2pos = d3.tsv.parseRows(pretext) ;
                        ms2pos = ms2pos.filter(function(d,unused_a,unused_b){return all_files.indexOf(d[0]) > -1 ;})
                        //Add bogus score as bogus field... terrible terrible hack...
                        ms2pos.forEach(function(d){d[9] = -1.0}) ;
                    }
                    complete_display() ;
                    d3.select("#animation").transition().duration(500).style('opacity',0.0) ;
                    d3.select("#animation").style("display","none") ;
                    d3.select("#launch_image").attr("src","img/dice_still.gif") ;
                    d3.select("#diceintitle").attr("src","img/dice_icon_small.png") ;
                }) ;
            }) ;
        } ;

        if(all_files){
            after_initialize_times() ;
        }else{
            d3.json("published",function(the_data){
                published = the_data ;
                active_category = published.variables[0].name ;
                all_files = [] ;
                published.files.forEach(function(pubfile){
                    published[pubfile.fname] = pubfile ;
                    all_files.push(pubfile.fname) ;
                }) ;
                published.variables.forEach(function(category){
                    categories[category.name] = category.values ;
                    category.values.forEach(function(val){
                        multifilter[val] = 1 ;
                    }) ;
                })
                d3.select("#filedomain").html(all_files.length) ;
                initialize_times(all_files,after_initialize_times) ;
            }) ;
        }
    };

    var score_timeout_var = setTimeout(score_all_ms2s,50) ;

    viewer = {
        "change_active_category": change_active_category,
        "flip_category": flip_category,
        "clear_variable": clear_variable,
        "clear_multifilter": clear_multifilter,
        "smz": smz,
        "load_em_all": load_em_all,
        "pset2prots": pset2prots,
        "protentry2prot": protentry2prot,
        "prot2peps": prot2peps,
        "proteo2peps": proteo2peps,
        "pep2mass": pep2mass,
        "pep2mass_non_event": pep2mass_non_event,
        "expose_ms2_frame": expose_ms2_frame,
        "expose_nist_frame": expose_nist_frame,
        "return_launch": return_launch
    }
    return viewer ;

}
