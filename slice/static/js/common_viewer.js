function initialize(color_by_source,file_filter,tolerance_by_source,shiftfilter,altfilter,labeler,metsets,metset_names){

    var metsets = metsets ;
    met_lookup = {}
    for(var i = 0 ; i < metset_names.length ; i++){
        var a_metset = metset_names[i] ;
        for(var j = 0; j < metsets[a_metset].length;j++){
            a_met = metsets[a_metset][j] ;
            met_lookup[a_met.name] = a_met ;
        }
    }
    console.log(met_lookup) ;
    load_metsets(metset_names) ;
    metset(metset_names[0]) ;

    function url_domain(data) {
        var    a      = document.createElement('a');
        a.href = data;
        return a.hostname;
    }

    var ms1target = null;
    var ms1time = null ;
    var shiftfilter_value = null ;
    var altfilter_value = null ;


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

    var hw = l10(1.0 + 10.0*ppm) ;

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
                .append("option").text(function(d) { return d[1] + " " + d[0]; }).property("value",function(d){return d[0]});
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


    function metset(the_metset){
        if(!the_metset){
            the_metset = d3.select("#metset").property("value")
        }
        console.log(the_metset) ;
        d3.select("#metabolite").selectAll("option").remove() ;
        d3.select("#metabolite").selectAll("option").data(metsets[the_metset]).enter()
                .append("option").text(function(d) { return d.name; }).property("value",function(d){return d.name});
    }






    function load_metsets(metlist){
        d3.select("#metset").selectAll("option").remove() ;
        d3.select("#metset").selectAll("option").data(metlist).enter()
                .append("option").text(function(d) { return d; }).property("value",function(d){return d ;});
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

    function met(){
        //var the_metset = d3.select("#metset").property("value") ;
        var the_met = d3.select("#metabolite").property("value") ;
        var details = met_lookup[the_met] ;
        d3.select("#t").property("value",details.rt) ;
        d3.select("#mz").property("value",details.mz) ;
        load_em_all() ;
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

    function xic(target_file,target_mz,callback){
        var the_vals = pic_data[target_file] ;
        var ts = file_times.get(target_file) ;
        answer = [] ;
        for(var i = 0; i < ts.length ; i++){
            answer.push({"time": ts[i], "intensity": the_vals[i][1]}) ;
        }
        xics.set(target_file,{"file": target_file, "mz": target_mz, "values": answer})
        callback() ;
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


//var color = d3.scale.category10();
    var color = color_by_source() ;

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

    //Write code to put stub svg in the ms1_viewer frame!!!
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

//target.append("line").attr("display", "none")
//        .attr("class", "targety")
//        .attr("x1", width - 6)
//        .attr("x2", width + 6);

    target.append("circle").attr("display", "none")
        .attr("class", "targetcircle").attr("r", 3.5) ;

    target.append("text").attr("id","targetlabel").text("").attr("display", "none");

    target.on("click",function(d){
        var altshift = false ;
        if (d3.event.shiftKey) {
            altshift = true ;
            if(shiftfilter_value){
                shiftfilter_value = null ;
            }else{
                if(ms1target){
                    console.log("ms1 -- filter by fraction") ;
                    console.log( shiftfilter(ms1target.file) ) ;
                    shiftfilter_value = shiftfilter(ms1target.file) ;
                }
            }
        }
        if (d3.event.altKey) {
            altshift = true ;
            if(altfilter_value){
                altfilter_value = null ;
            }else{
                if(ms1target){
                    console.log("filter by patient") ;
                    console.log( altfilter(ms1target.file) ) ;
                    altfilter_value = altfilter(ms1target.file) ;
                }
            }
        }
        if(altshift){
            transdisplay()
        }else{
            // target_domain = url_domain(document.URL) ;
            // console.log(target_domain) ;
            // prefix = target_domain.slice(0,4) ;
            // if(prefix == "beta"){
            //     target_domain = "betaslice" + target_domain.slice(8) ;
            // }else{
            //     target_domain = "slice" + target_domain.slice(4) ;
            // }
            // dataset = document.URL.split("/").slice(-2,-1)[0] ;
            // https = document.URL.slice(0, document.URL.indexOf(":"))
            // window.open(https + "://" + target_domain + "/" + dataset + "/files/" + ms1target.file +  "/scans/" + ms1time) ;
            window.open(document.URL.slice(0,document.URL.search("xic"))+"files/" + ms1target.file +  "/scans/" + ms1time) ; 
        }


//        if(ms1target){
//            render_precursor(d3.select("#mz").property("value"),ms1target.file,ms1time) ;
//            d3.select('#ms1_holder').style('display','inline') ;  console.log(ms1target.file)
//        }
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
        //target.select(".targety").attr("display", "none");
    }

    function mousemove() {

        var in_x = d3.mouse(this)[0] ;
        var in_y = d3.mouse(this)[1] ;

        var the_x = x.invert(d3.mouse(this)[0]) ;
        var the_y = y.invert(d3.mouse(this)[1]) ;

        var bisect = d3.bisector(function(d) { return d.time; }).right;
        //var d = offsets[Math.round((x.invert(d3.mouse(this)[0]) - start) / step)];
        var candidates = [] ;
        xics.forEach(function(fname,xc){
            if(shiftfilter_value){
                if(shiftfilter(fname) != shiftfilter_value){
                    return ;
                }
            }
            if(altfilter_value){
                if(altfilter(fname) != altfilter_value){
                    return ;
                }
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
            //target.select(".targety").attr("transform", "translate(0," + c_y + ")");
            target.select(".targetcircle").attr("display", "inline") ;
            target.select("#targetlabel").text(labeler(closest.file)).attr("fill","black").attr("transform", "translate(" + (c_x + 10) + "," + (c_y +5) + ")").attr("display", "inline") ;;
            target.select(".targetx").attr("display", "inline");
        } else {
            ms1target = null ;
            target.select(".targetcircle").attr("display", "none");
            target.select("#targetlabel").attr("display","none") ;
            target.select(".targetx").attr("display", "none");
        }
        //target.select(".targety").attr("display", "inline");

        //svg.selectAll(".x.axis path").style("fill-opacity", Math.random()); // XXX Chrome redraw bug
    }

    var score_threshold = 0.2 ;

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
//    x.domain([
//        d3.min(xics.values(), function(an_xic){ return d3.min(an_xic.values, function(v){ return v.time; });}) ,
//        d3.max(xics.values(), function(an_xic){ return d3.max(an_xic.values, function(v){ return v.time; });})
//    ]) ;
//
//    y.domain([
//        d3.min(xics.values(), function(an_xic){ return d3.min(an_xic.values, function(v){ return v.intensity; });}) ,
//        d3.max(xics.values(), function(an_xic){ return d3.max(an_xic.values, function(v){ return v.intensity; });})
//    ]) ;

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

        //var brush_start = y_max_time - 8 ;
        //var brush_stop = y_max_time + 8 ;

        var brush_start = target_time_start ;
        var brush_stop = target_time_stop ;

        x.domain([brush_start,brush_stop]) ;
        y.domain([0,y_max]) ;

        var area2 = d3.svg.area().interpolate("monotone").x(function(d){return x2(d.time)}).y0(height2).y1(function(d){return y2(d.intensity)}) ;
        context.append("path").attr("id","bpc").attr("d",area2(the_bpc)) ;

        var xic_elems = focus.selectAll(".xic").data(xics.values()) ;
        xic_elems.enter().append("g").attr("class","xic").append("path").attr("class","line").attr("clip-path", "url(#clip)") ;
//    xic_elems.select("path").attr("d", function(d){return line(d.values);})
//            .style("stroke", function(d){return color(d.file);})
//            ;

        function score_blurb(a_score){
            if(a_score < 0){
                return "" ;
            }else{
                return " (" + a_score.toFixed(2) + ")" ;
            }
        }

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

        var grand_maxer = function(v){
            if((v.time > x.domain()[0]) && (v.time < x.domain()[1])){
            }
        }

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
            if(shiftfilter_value){
                if(shiftfilter(d.file)!=shiftfilter_value){
                    return ;
                }
            }
            if(altfilter_value){
                if(altfilter(d.file)!=altfilter_value){
                    return ;
                }
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
//        return line(d.values);
            if(shiftfilter_value){
                if(shiftfilter(d.file)!=shiftfilter_value){
                    //Should be the same as below...
                    return line([]);
                    //return line([{"time": x.domain()[0], "intensity": 0},{"time": x.domain()[1], "intensity": 0}]);
                }
            }
            if(altfilter_value){
                if(altfilter(d.file)!=altfilter_value){
                    //Should be the same as below...
                    return line([]) ;
                    //return line([{"time": x.domain()[0], "intensity": 0},{"time": x.domain()[1], "intensity": 0}]);
                }
            }
            var subdata = dfilter(d.values) ;
            if(subdata.length == 0){
                //Should be the same as above...
                return line([]) ;
                //return line([{"time": x.domain()[0], "intensity": 0},{"time": x.domain()[1], "intensity": 0}]);
            } else {
                return line(subdata);
            }
        }).style("stroke", function(d){return color(d.file);});

        target.selectAll(".point")
            .attr("x", function(d, i) { return x(d[1]) })
            .attr("y", function(d, i) { return "-10px" })
            .attr("display", function(d, i){

            if(shiftfilter_value){
                if(shiftfilter(d[0])!=shiftfilter_value){
                    //Should be the same as below...
                    return "none" ;
                }
            } ;

            if(altfilter_value){
                if(altfilter(d[0])!=altfilter_value){
                    //Should be the same as below...
                    return "none" ;
                }
            } ;


                if( (d[1] > x.domain()[0]) && (d[1] < x.domain()[1])){return "inline"}else{return "none"}}
            ) ;

        focus.select(".x.axis").call(xAxis) ;
        focus.select(".y.axis").call(yAxis) ;

        hidetarget() ;

        var shifters = d3.set() ;
        var alters = d3.set() ;
        visibles.forEach(function(vis){
            shifters.add(shiftfilter(vis)) ;
            alters.add(altfilter(vis)) ;
        }) ;
        d3.select("#totals").html(shiftfilter.label + "=" + shifters.values().length + " " + altfilter.label + "=" + alters.values().length + " Files=" + visibles.values().length) ;
        inject_report() ;
    }


    function transdisplay(){

    var grand_max = 0 ;
    var visibles = d3.set() ;

    var grand_maxer = function(v){
        if((v.time > x.domain()[0]) && (v.time < x.domain()[1])){
        }
    }

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
        if(shiftfilter_value){
            if(shiftfilter(d.file)!=shiftfilter_value){
                return ;
            }
        }
        if(altfilter_value){
            if(altfilter(d.file)!=altfilter_value){
                return ;
            }
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
//        return line(d.values);
        if(shiftfilter_value){
            if(shiftfilter(d.file)!=shiftfilter_value){
                //Should be the same as below...
                return line([]);
                //return line([{"time": x.domain()[0], "intensity": 0},{"time": x.domain()[1], "intensity": 0}]);
            }
        }
        if(altfilter_value){
            if(altfilter(d.file)!=altfilter_value){
                //Should be the same as below...
                return line([]) ;
                //return line([{"time": x.domain()[0], "intensity": 0},{"time": x.domain()[1], "intensity": 0}]);
            }
        }
        var subdata = dfilter(d.values) ;
        if(subdata.length == 0){
            //Should be the same as above...
            return line([]) ;
            //return line([{"time": x.domain()[0], "intensity": 0},{"time": x.domain()[1], "intensity": 0}]);
        } else {
            return line(subdata);
        }
    }).style("stroke", function(d){return color(d.file);});

    target.selectAll(".point")
        .attr("x", function(d, i) { return x(d[1]) })
        .attr("y", function(d, i) { return "-10px" })
        .attr("display", function(d, i){

        if(shiftfilter_value){
            if(shiftfilter(d[0])!=shiftfilter_value){
                //Should be the same as below...
                return "none" ;
            }
        } ;

        if(altfilter_value){
            if(altfilter(d[0])!=altfilter_value){
                //Should be the same as below...
                return "none" ;
            }
        } ;


            if( (d[1] > x.domain()[0]) && (d[1] < x.domain()[1])){return "inline"}else{return "none"}}
        ) ;

    focus.select(".x.axis").call(xAxis) ;
    focus.select(".y.axis").call(yAxis) ;

    hidetarget() ;

    var shifters = d3.set() ;
    var alters = d3.set() ;
    visibles.forEach(function(vis){
        shifters.add(shiftfilter(vis)) ;
        alters.add(altfilter(vis)) ;
    }) ;
    d3.select("#totals").html(shiftfilter.label + "=" + shifters.values().length + " " + altfilter.label + "=" + alters.values().length + " Files=" + visibles.values().length) ;
    inject_report() ;
}



    function inject_report(){
        xic_map = d3.map() ;
        scan_map = d3.map() ;

        if(xics.keys().length < all_files.length){
            return ;
        }

        xics.keys().forEach(function(k){
            var catkey = color.label(k) ;
            if(!xic_map.has(catkey)){
                xic_map.set(catkey,0)
            }
            xic_map.set(catkey,xic_map.get(catkey)+1)
        }) ;

        report_html = "" ;
        xic_map.keys().forEach(function(k){
            if(scan_map.get(k)){
                if(d3.select("#offscore").property("checked") == true){
                report_html =
                    report_html +
                        "<tr><td style='text-align: right;font-weight: bold;font-style: sans-serif; color: " + color.label_color(k) + "'>" + k + "</td><td>" +
                        ": " + xic_map.get(k) + " xics, " + scan_map.get(k)[0] + " ms2 scans</td></tr>"
                }else {
                report_html =
                    report_html +
                        "<tr><td style='text-align: right;font-weight: bold;font-style: sans-serif; color: " + color.label_color(k) + "'>" + k + "</td><td>" +
                        ": " + xic_map.get(k) + " xics, " + scan_map.get(k)[0] + " ms2 scans, " + scan_map.get(k)[1] + " scored, yielding " + scan_map.get(k)[2] + " hits</td></tr>"

                }
            }else{
                report_html =
                    report_html +
                        "<tr><td style='text-align: right;font-weight: bold;font-style: sans-serif; color: " + color.label_color(k) + "'>" + k + "</td><td>" +
                        ": " + xic_map.get(k) + " xics, m/z=" + hidden_mz + ", time=" + ((target_time_start+target_time_stop)/2.0).toFixed(2) +"</td></tr>"
            }
        }) ;
        d3.select("#report_table").html(
            report_html
//            "<tr><td style='text-align: right;font-weight: bold;color: purple'>P5</td><td>"
//            + ": 160 xics, 210 ms2 scans, 10 scored, yielding 3 hits</td></tr>"
//            + "<tr><td style='text-align: right;font-weight: bold;color: orange'>P6</td><td>"
//            + ": 150 xics, 210 ms2 scans, 10 scored, yielding 3 hits</td></tr>"
        )
    }

    function animLoop( render ) {
        var running, lastFrame = +new Date,
            raf = window.mozRequestAnimationFrame    ||
                window.webkitRequestAnimationFrame ||
                window.msRequestAnimationFrame     ||
                window.oRequestAnimationFrame;
        function loop( now ) {
            // stop the loop if render returned false
            if ( running !== false ) {
                raf( loop );
                var deltaT = now - lastFrame;
                //running = render()
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
        color.domain(file_list);
        file_list.forEach(function(fname){
            var the_times = [] ;
            var the_xic = pic_data[fname] ;
            for(var i = 0 ; i < the_xic.length ; i++){
                the_times[i] = the_xic[i][0] ;
            }
            file_times.set(fname,the_times) ;
        });
        finisher() ;
    }


    var all_files = undefined ;

    var ms2pos = undefined ;
    var hidden_mz = undefined ;
    var hidden_seq = undefined ;
    var hidden_charge = undefined ;

    var pic_data = null ;
    var target_time_start = null ;
    var target_time_stop = null ;

    function load_em_all(mz_query,all_data,t_start,t_stop){
        if(!mz_query){
            mz_query = parseFloat(d3.select("#mz").property("value")).toFixed(4) ;
            var t_query_val = parseFloat(d3.select("#t").property("value")) ;
            var t_query = t_query_val.toFixed(2) ;
            tq_start = t_query_val-0.5 ;
            tq_stop = t_query_val+0.5 ;
            d3.select("#launch_image").attr("src","../../../img/dice_rotate.gif") ;
            d3.json(document.URL.slice(0,document.URL.search("xic"))+"pic/"+mz_query+"/"+ t_query, function(error, json) {
                if (error) return console.warn(error);
                load_em_all(mz_query,json,tq_start,tq_stop);
            });
            return ;
        }
        hidden_mz = mz_query ;
        hidden_seq = d3.select("#metabolite").property("value") ;


        if(!all_data){
            return ;
        }

        console.log(all_data) ;
        pic_data = all_data["ric_data"] ;
        file_list = all_data["file_list"] ;
        console.log(file_list) ;
        target_time_start = t_start ;
        target_time_stop = t_stop ;


        d3.select("#launch_image").attr("src","../../../img/dice_rotate.gif") ;

        var after_initialize_times = function(){
            xics = d3.map() ;
            all_files.forEach(function(fname){
                xic(fname,mz_query,function(){
                    if(xics.keys().length == all_files.length){
                        complete_display() ;
                        d3.select("#launch_image").attr("src","../../../img/dice_still.gif") ;
                    }
                })
            }) ;
        } ;

        if(all_files){
            after_initialize_times() ;
        }else{
                all_files = file_list ;
                all_files = all_files.filter(file_filter) ;
                d3.select("#filedomain").html(all_files.length) ;
                all_files.sort(function(a,b){
                    var label_a = color.label(a) ;
                    var label_b = color.label(b) ;
                    return color.sort_order.indexOf(label_a) - color.sort_order.indexOf(label_b) ;
                })
                initialize_times(all_files,after_initialize_times) ;
        }
    };

    viewer = {
        "smz": smz,
        "load_em_all": load_em_all,
        "metset": metset,
        //"protentry2prot": protentry2prot,
        //"prot2peps": prot2peps,
        //"proteo2peps": proteo2peps,
        "met": met,
        //"pep2mass_non_event": pep2mass_non_event,
        "return_launch": return_launch
    }
    return viewer ;

}
