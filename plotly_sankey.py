from plotly import graph_objects
import sys
import math

MAX=20000000

def to_rgba(color):
    red = int(color[1:3], 16)
    green = int(color[3:5], 16)
    blue = int(color[5:7], 16)
    alpha = int(color[7:9], 16)/255.0
    return f"rgba({red}, {green}, {blue}, {alpha})"

def threshold(num):
    return num if num<MAX else MAX

def sankey(dag, node_properties, title):
    # reformat DAG into inputs for Plotly
    labels_to_indices = {}
    labels = []
    node_colors = []
    sources = []
    targets = []
    values = []
    colors = []
    x = []
    y = []
    

    for node in dag:
        if node not in labels:
            labels.append(node)
            labels_to_indices[node] = len(labels)-1
            node_colors.append(node_properties[node]["color"])
        for target_node in dag[node]:
            if target_node not in labels:
                labels.append(target_node)
                labels_to_indices[target_node] = len(labels)-1
                node_colors.append(node_properties[target_node]["color"])
            sources.append(labels_to_indices[node])
            targets.append(labels_to_indices[target_node])
            values.append(dag[node][target_node]['weight'])
            colors.append(dag[node][target_node]['color'])
    
    # create Sankey plot
    fig = graph_objects.Figure(data=[
        graph_objects.Sankey(
            arrangement="fixed",
            node = {
                "pad": 15,
                "thickness": 35,
                "line": {"color": "black", "width": 0.5},
                "label": labels,
                "color": node_colors
            },
            link = {
                "source": sources,
                "target": targets,
                "value": values,
                "color": colors
            },
            textfont = {
                "size": 20
            }
            )
        ]
    )
    fig.update_layout(title_text=title, font_size=10)

    return fig

file_data = {}

with open(sys.argv[1]) as f:
    for line in f:
        rec = line.split("\t")
        if rec[1] == "pct_miscalls" or rec[1] == "pct_clinvar_miscalls":
            continue
        file_data[".".join(rec[1:4]).strip()] = int(rec[0])

dag_dict_all = {
    "All Cs": 
    {
        "High Mappability": {"weight": (file_data["calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.nofilt"]), "color": to_rgba("#2A2A6680")}, 
        "Low Mappability": {"weight": (file_data["lowmap_calls.CpG.nofilt"]), "color": to_rgba("#80808080")}, 
        "Low MapQ": {"weight": (file_data["calls_mq0.CpG.nofilt"]-file_data["calls.CpG.nofilt"]), "color": to_rgba("#BA5B2C80")}
    }, 
    "High Mappability": 
    {
        "Misfiltered": {"weight": (file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#99553380")}, 
        "Callable Cs": {"weight": ((file_data["calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.nofilt"])-file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#3333CC80")}
    },
    "Low Mappability": 
    {
        "Rescued": {"weight": (file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#33339980")}, 
        "Removed": {"weight": (file_data["lowmap_calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#CC663380")}
    },
    "Low MapQ": 
    {
        "Uncallable Cs": {"weight": (file_data["calls_mq0.CpG.nofilt"]-file_data["calls.CpG.nofilt"]), "color": to_rgba("#FF7D3C80")}
    },
    "Removed": 
    {
        "Uncallable Cs": {"weight": (file_data["lowmap_calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#FF7D3C80")}
    },
    "Misfiltered": 
    {
        "Uncallable Cs": {"weight": (file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#FF7D3C80")}
    },
    "Rescued": 
    {
        "Callable Cs": {"weight": (file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#3333CC80")}
    }
}


dag_dict_focused = {
    "Affected Cs": 
    {
        "High Mappability": {"weight": (file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#2A2A6680")}, 
        "Low Mappability": {"weight": (file_data["lowmap_calls.CpG.nofilt"]), "color": to_rgba("#80808080")} 
    }, 
    "High Mappability": 
    {
        "Misfiltered": {"weight": (file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#99553380")} 
    },
    "Low Mappability": 
    {
        "Rescued": {"weight": (file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#33339980")}, 
        "Removed": {"weight": (file_data["lowmap_calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#CC663380")}
    },
    "Removed": 
    {
        "Uncallable Cs": {"weight": (file_data["lowmap_calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#FF7D3C80")}
    },
    "Misfiltered": 
    {
        "Uncallable Cs": {"weight": (file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#FF7D3C80")}
    },
    "Rescued": 
    {
        "Callable Cs": {"weight": (file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#3333CC80")}
    }
}

dag_dict_thresholded = {
    "All Cs": 
    {
        "High Mappability": {"weight": threshold(file_data["calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.nofilt"]), "color": to_rgba("#2A2A6680")}, 
        "Low Mappability": {"weight": threshold(file_data["lowmap_calls.CpG.nofilt"]), "color": to_rgba("#80808080")}, 
        "Low MapQ": {"weight": threshold(file_data["calls_mq0.CpG.nofilt"]-file_data["calls.CpG.nofilt"]), "color": to_rgba("#BA5B2C80")}
    }, 
    "High Mappability": 
    {
        "Misfiltered": {"weight": threshold(file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#99553380")}, 
        "Callable Cs": {"weight": threshold((file_data["calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.nofilt"])-file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#3333CC80")}
    },
    "Low Mappability": 
    {
        "Rescued": {"weight": threshold(file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#33339980")}, 
        "Removed": {"weight": threshold(file_data["lowmap_calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#CC663380")}
    },
    "Low MapQ": 
    {
        "Uncallable Cs": {"weight": threshold(file_data["calls_mq0.CpG.nofilt"]-file_data["calls.CpG.nofilt"]), "color": to_rgba("#FF7D3C80")}
    },
    "Removed": 
    {
        "Uncallable Cs": {"weight": threshold(file_data["lowmap_calls.CpG.nofilt"]-file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#FF7D3C80")}
    },
    "Misfiltered": 
    {
        "Uncallable Cs": {"weight": threshold(file_data["underfiltered_calls.CpG.combined"]), "color": to_rgba("#FF7D3C80")}
    },
    "Rescued": 
    {
        "Callable Cs": {"weight": threshold(file_data["lowmap_calls.CpG.filt"]), "color": to_rgba("#3333CC80")}
    }
}

node_properties = {
    "All Cs": {"color": "#808080"}, # "x": 0, "y": 0},
    "Affected Cs": {"color": "#808080"}, # "x": 0, "y": 0},
    "High Mappability": {"color": "#2A2A66"}, # "x": 0.2, "y": 0}, 
    "Low Mappability": {"color": "#808080"}, # "x": 0.2, "y": 0.733}, 
    "Low MapQ": {"color": "#BA5B2C"}, # "x": 0.6, "y": 1.0}, 
    "Misfiltered": {"color": "#995533"}, # "x": 0.8, "y": 0.8}, 
    "Callable Cs": {"color": "#3333CC"}, # "x": 1.0, "y": 0}, 
    "Rescued": {"color": "#333399"}, # "x": 0.5, "y": 0.8}, 
    "Removed": {"color": "#CC6633"}, # "x": 0.5, "y": 0.4}, 
    "Uncallable Cs": {"color": "#FF7D3C"}, # "x": 1.0, "y": 0.7}, 
    }
fig1 = sankey(dag_dict_all, node_properties, "C Filtering Overview")
fig2 = sankey(dag_dict_focused, node_properties, "C Filtering Overview")
fig3 = sankey(dag_dict_thresholded, node_properties, "C Filtering Overview")

fig1.write_image("sankey_all_"+sys.argv[2]+".svg", width=1200, height=900)
fig2.write_image("sankey_focused_"+sys.argv[2]+".svg", width=1200, height=900)
fig3.write_image("sankey_thresholded_"+sys.argv[2]+".svg", width=1200, height=900)
