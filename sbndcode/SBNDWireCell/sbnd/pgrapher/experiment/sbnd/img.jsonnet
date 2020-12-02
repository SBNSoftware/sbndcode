// some functions to help build pipelines for imaging.  These are
// mostly per-apa but tiling portions are per-face.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

{
    // A functio that sets up slicing for an APA.
    slicing :: function(anode, aname, tag="", span=4) {
        ret: g.pnode({
            type: "SumSlices",
            name: "slicing-"+aname,
            data: {
                tag: tag,
                tick_span: span,
                anode: wc.tn(anode),
            },
        }, nin=1, nout=1, uses=[anode]),
    }.ret,

    // A function sets up tiling for an APA incuding a per-face split.
    tiling :: function(anode, aname) {

        local slice_fanout = g.pnode({
            type: "SliceFanout",
            name: "slicefanout-" + aname,
            data: { multiplicity: 2 },
        }, nin=1, nout=2),

        local tilings = [g.pnode({
            type: "GridTiling",
            name: "tiling-%s-face%d"%[aname, face],
            data: {
                anode: wc.tn(anode),
                face: face,
            }
        }, nin=1, nout=1, uses=[anode]) for face in [0,1]],

        local blobsync = g.pnode({
            type: "BlobSetSync",
            name: "blobsetsync-" + aname,
            data: { multiplicity: 2 }
        }, nin=2, nout=1),

        ret: g.intern(
            innodes=[slice_fanout],
            outnodes=[blobsync],
            centernodes=tilings,
            edges=
                [g.edge(slice_fanout, tilings[n], n, 0) for n in [0,1]] +
                [g.edge(tilings[n], blobsync, 0, n) for n in [0,1]],
            name='tiling-' + aname),
    }.ret,

    // Just clustering
    clustering :: function(anode, aname, spans=1.0) {
        ret : g.pnode({
            type: "BlobClustering",
            name: "blobclustering-" + aname,
            data:  { spans : spans }
        }, nin=1, nout=1),
    }.ret, 

    // this bundles clustering, grouping and solving.  Other patterns
    // should be explored.  Note, anode isn't really needed, we just
    // use it for its ident and to keep similar calling pattern to
    // above..
    solving :: function(anode, aname, spans=1.0, threshold=0.0) {
        local bc = g.pnode({
            type: "BlobClustering",
            name: "blobclustering-" + aname,
            data:  { spans : spans }
        }, nin=1, nout=1),
        local bg = g.pnode({
            type: "BlobGrouping",
            name: "blobgrouping-" + aname,
            data:  {
            }
        }, nin=1, nout=1),
        local bs = g.pnode({
            type: "BlobSolving",
            name: "blobsolving-" + aname,
            data:  { threshold: threshold }
        }, nin=1, nout=1),
        ret: g.intern(
            innodes=[bc], outnodes=[bs], centernodes=[bg],
            edges=[g.edge(bc,bg), g.edge(bg,bs)],
            name="solving-" + aname),
    }.ret,

    dump :: function(anode, aname, drift_speed) {
        local js = g.pnode({
            type: "JsonClusterTap",
            name: "clustertap-" + aname,
            data: {
                filename: "clusters-"+aname+"-%04d.json",
                drift_speed: drift_speed
            },
        }, nin=1, nout=1),

        local cs = g.pnode({
            type: "ClusterSink",
            name: "clustersink-"+aname,
            data: {
                filename: "clusters-apa-"+aname+"-%d.dot",
            }
        }, nin=1, nout=0),
        ret: g.intern(innodes=[js], outnodes=[cs], edges=[g.edge(js,cs)],
                      name="clusterdump-"+aname)
    }.ret,

    // A function that reverts blobs to frames
    reframing :: function(anode, aname) {
        ret : g.pnode({
            type: "BlobReframer",
            name: "blobreframing-" + aname,
            data: {
                frame_tag: "reframe%d" %anode.data.ident,
            }
        }, nin=1, nout=1),
    }.ret,

    // fill ROOT histograms with frames
    magnify :: function(anode, aname, frame_tag="orig") {
        ret: g.pnode({
          type: 'MagnifySink',
          name: 'magnify-'+aname,
          data: {
            output_filename: "magnify-img.root",
            root_file_mode: 'UPDATE',
            frames: [frame_tag + anode.data.ident],
            trace_has_tag: true,
            anode: wc.tn(anode),
          },
        }, nin=1, nout=1),
    }.ret,

    // the end
    dumpframes :: function(anode, aname) {
        ret: g.pnode({
            type: "DumpFrames",
            name: "dumpframes-"+aname,
        }, nin=1, nout=0),
    }.ret,

}
