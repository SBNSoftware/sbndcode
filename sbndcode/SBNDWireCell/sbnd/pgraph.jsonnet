// This file provides various helper functions that assist in
// configuring jobs which take a graph of nodes, eg Pgrapher.  The
// nomenclature here is that an "inode" is a configuration object
// corresponding to an WCT INode component and a "pnode" is a wrapper
// around an inode or a number of pnodes which assist in constructing
// the graph in a piecewise manner.  It is best to only make pnodes
// via function calls, otherwise some functionality may not work.

local wc = import "wirecell.jsonnet";

{
    // Construct a port structure as used to form 1/2 of an edge.
    port(inode, num=0) :: { node: wc.tn(inode), port: num },


    // Make an edge between two pnodes by passing those pnodes as objects
    edge(tail, head, tp=0, hp=0):: {
        tail: tail.oports[tp],
        head: head.iports[hp],
    },

    // make an edge by passing two pnode "type:name" labels and
    // optional port numbers.
    edge_labels(tlabel, hlabel, tp=0, hp=0):: {
        tail: {
            node: tlabel,
            port: tp,
        }, 
        head: {
            node: hlabel,
            port: hp,
        },
    },
    

    // Break an existing edge, terminating the tail end with a new
    // head and starting the head end with a new tail.  Graphically:
    // edge --> [edge[tail]->nh, nt->edge[head]]
    break_edge(edge, nh, nt):: [
        {
            tail: edge.tail,
            head: nh
        },
        {
            tail: nt,
            head: edge.head,
        },
    ],


    // Break and existing edge at the index in the edges array, return
    // new array of edges with the inserted new head and tail ports.
    break_insert_edge(index, edges, nh, nt)::
    std.join($.break_edge(edges[index], nh, nt), [edges[0:index], edges[index+1:std.length(edges)]]),


    // Strip any pnodes    
    strip_pnodes(arr):: std.filter(function(x) x.type != "Pnode", arr),

    // Return true if is not null/empty
    isSomething(val) ::
    if val == null then false
    else if std.type(val) == "array" then std.length(val) > 0
    else if std.type(val) == "object" then std.length(val) > 0
    else true,

    // Return a new object with key removed, if it exists.
    // Note, there is std.prune() but it's pretty slow on big objects as it recurs
    prune_key(obj, key) ::
    if std.objectHas(obj, key)
    then {
        [k]: obj[k]
        for k in std.objectFields(obj) if k != key
    }
    else obj,

    // Return a new list made from the input with any elmenets which are null or empty lists removed
    prune_array(arr) :: [ x for x in arr if $.isSomething(x) ],

    // Helper recursively find all objects in "uses" array, removing
    // the array asit goes.  Return catenation of list "l" and all "uses" found.
    popuses(l, obj):: if std.objectHas(obj, 'uses')
    then l + std.foldl($.popuses, obj.uses, []) + [$.prune_key(obj, 'uses')]
    else l + [obj],

    // Return all "uses" objects.  Note, the returned list may need to
    // be passed to wc.unique_list().
    resolve_uses(seq):: $.strip_pnodes(std.foldl($.popuses, seq, [])),
    
    // Make a pnode from an inode, giving its input and output port
    // multiplicity.  Use this instead of creating a pnode by hand.
    // Any other WCT component objects which are referenced by this
    // one should be passed in "uses" (or, as a special inode.uses).
    pnode(inode, nin=0, nout=0, uses=[]):: {
        type: "Pnode",
        name: wc.tn(inode),
        edges: [],
        uses: uses + [inode],
        iports: [$.port(inode, n) for n in std.range(0,nin)][:nin],
        oports: [$.port(inode, n) for n in std.range(0,nout)][:nout],
    },


    // Produce a new pnode from collections of input and output pnodes
    // and any internal nodes and edges.  The resulting "uses" and
    // "edges" are then resolved, aggregated, flattened.  Unless
    // explicitly given, all iports of innodes become iports of the
    // new pnode, etc for output.
    intern(innodes=[], outnodes=[], centernodes=[], edges=[], iports=[], oports=[], name=""):: {
        local nodes = innodes+outnodes+centernodes,
        type: "Pnode",
        name: name,
        uses: nodes,
        edges: $.prune_array(edges + std.flattenArrays([n.edges for n in nodes])),
        iports: if std.length(iports) == 0 then std.flattenArrays([n.iports for n in innodes]) else iports,
        oports: if std.length(oports) == 0 then std.flattenArrays([n.oports for n in outnodes]) else oports,
    },


    // Intern an ordered list of elements to form a linear pipeline by
    // subsequently connecting one element's output port 0 to the next
    // element's input port 0.  The first/last elements iport/oport
    // will be used for the pipeline's iport/oport if existing.
    pipeline(elements, name=""):: {
        local nele = std.length(elements),
        local pedges = [$.edge(elements[i], elements[i+1]) for i in std.range(0,nele-2)],
        type: "Pnode",
        name: name,
        uses: elements,
        edges: $.prune_array(pedges + std.flattenArrays([n.edges for n in elements])),
        iports: if std.length(elements[0].iports) == 0 then [] else [elements[0].iports[0]],
        oports: if std.length(elements[nele-1].oports) == 0 then [] else [elements[nele-1].oports[0]],
    },


    // Return a new pnode built by breaking an existing edge at given
    // index and patching the break with the given head and tail nodes
    // and their ports.  If a name is given it is set else the name of
    // the original pnode is kept.
    insert_one(pnode, index, newhead, newtail, iport=0, oport=0, name=null):: {
        type: "Pnode",
        name: $.prune_array([name, pnode.name])[0],
        uses: [pnode,newhead,newtail],
        edges: $.break_insert_edge(index, pnode.edges, newhead.iports[iport], newtail.oports[oport]) + newhead.edges + newtail.edges,
        iports: pnode.iports,
        oports: pnode.oports,
    },

    // Return a list of indices where item is found in list
    find_indices(list, item):: std.filter(std.isNumber, 
                                           std.mapWithIndex(function(ind,ele)
                                                            if ele == item
                                                            then ind
                                                            else null,
                                                            list)),

    // Like insert_one() but give edge to break instead of index
    insert_node(pnode, edge_to_break, newhead, newtail, iport=0, oport=0, name=null):: 
    self.insert_one(pnode, self.find_indices(pnode.edges, edge_to_break)[0], newhead, newtail, iport, oport, name),


    // Joint N sources using joiner, return pnode that looks like a
    // single source.  The joiner must be capable of handling and
    // N-join.  Each source is connected to joiner's input ports in
    // order.
    join_sources(joiner, sources, n=2) :: $.intern(outnodes=[joiner],
                                                   centernodes=sources,
                                                   iports=[],
                                                   edges=std.mapWithIndex(function(ind,s) $.edge(s,joiner,0,ind),
                                                                          sources),
                                                  ),
    

    // Call this to return the edges from a graph (a pnode).  It takes
    // care to remove any duplicates which can be slow so do NOT call
    // this except when getting a final list of edges.
    edges(graph) :: wc.unique_list(graph.edges),

    // Call this to return the final "uses" list which can be used as
    // part of the final wire cell configuration sequence.  It
    // recursively finds the uses of all uses (dawg) and returns a
    // unique list.  Do NOT call this except at high level as it's
    // somewhat expensive and need not be called on intermediate uses
    // lists.
    uses(graph) :: wc.unique_list(self.resolve_uses(graph.uses)),
    

    // Some utility functions to build fan-out/in subgraphs.
    fan:: {
        // Build a fanout-[pipelines]-fanin graph.  pipelines is a list of
        // pnode objects, one for each spine of the fan.
        pipe :: function(fout, pipelines, fin, name="fanpipe", outtags=[], tag_rules=[]) {

            local fanmult = std.length(pipelines),

            local fanout = $.pnode({
                type: fout,
                name: name,
                data: {
                    multiplicity: fanmult,
                    tag_rules: tag_rules,
                },
            }, nin=1, nout=fanmult),


            local fanin = $.pnode({
                type: fin,
                name: name,
                data: {
                    multiplicity: fanmult,
                    tags: outtags,
                },
            }, nin=fanmult, nout=1),

            ret: $.intern(innodes=[fanout],
                          outnodes=[fanin],
                          centernodes=pipelines,
                          edges=
                          [$.edge(fanout, pipelines[n], n, 0) for n in std.range(0, fanmult-1)] +
                          [$.edge(pipelines[n], fanin, 0, n) for n in std.range(0, fanmult-1)],
                          name=name),
        }.ret,

        // Build a fanout-[pipelines] graph where each pipe is self
        // terminated.  pipelines is a list of pnode objects, one for each
        // spine of the fan.
        sink :: function(fout, pipelines, name="fansink", tag_rules=[]) {

            local fanmult = std.length(pipelines),

            local fanout = $.pnode({
                type: fout,
                name: name,
                data: {
                    multiplicity: fanmult,
                    tag_rules: tag_rules,
                },
            }, nin=1, nout=fanmult),


            ret: $.intern(innodes=[fanout],
                          outnodes=[],
                          centernodes=pipelines,
                          edges=
                          [$.edge(fanout, pipelines[n], n, 0) for n in std.range(0, fanmult-1)],
                          name=name),
        }.ret,
    },                          // fan
}
