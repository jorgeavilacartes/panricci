def get_sources_sinks(path_gfa):
    sources = []
    sinks   = []
    with open(path_gfa, "r") as fp:
        for line in fp.readlines():

            # paths
            nodes_path=[]
            if line.startswith("P"):
                _, seq_id, path, *_ = line.replace("\n","").split("\t")

                nodes = path.split(",")
                source = nodes[0]
                sink = nodes[-1]

                sourceid = source.replace("+","").replace("-","")
                sinkid = sink.replace("+","").replace("-","")

                sources.append(sourceid)
                sinks.append(sinkid)
    
    sources, sinks = list(set(sources)), list(set(sinks))
    return sources, sinks
