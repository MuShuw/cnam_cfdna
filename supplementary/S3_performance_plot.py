from functools import reduce

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import seaborn as sn

from bokeh.layouts import column, row
from bokeh.models import Select, ColumnDataSource, HoverTool
from bokeh.palettes import Category20
from bokeh.plotting import curdoc, figure
from bokeh.transform import factor_cmap, factor_mark
from bokeh.models.widgets import CheckboxGroup
from bokeh.models.widgets.markups import Div

dt = pd.read_csv("S2_grouped_performance.tsv", sep="\t", index_col=0)
dt_acu=dt[dt["performance"]=="accuracy_mean"]

DT={}

METRICS=['mean','median','min','max']

for i in dt_acu["ontologyGroup"].unique():
    DT[i]=dt_acu[dt_acu["ontologyGroup"]==i]
    rename_dict={}
    for met in METRICS:
        rename_dict[met]='_'.join([i,met])
    DT[i]=DT[i].rename(rename_dict,axis=1).drop(["ontologyGroup","performance"], axis=1)

data_frames=list(DT.values())

DF=reduce(lambda left,right: pd.merge(left,right), data_frames)

df = DF.copy()

RNAS=list(df["rna_type"].unique())
GROUPS=list(df["group"].unique())
FEATS=list(df["features"].unique())
SORTS=list(df["sortingMethod"].unique())

SIZES = list(range(9, 40, 2))
COLORS = Category20[20]
N_SIZES = len(SIZES)
N_COLORS = len(COLORS)

columns = sorted(df.columns)
variable = [x for x in columns if any( i in x for i in METRICS)]
discrete = [x for x in columns if df[x].dtype == object]
continuous = [x for x in columns if x not in discrete]

def create_figure():
    df = DF[(DF["rna_type"].isin([RNAS[xx] for xx in rna_box.active])) &
            (DF["group"].isin([GROUPS[xx] for xx in group_box.active])) &
            (DF["features"].isin([FEATS[xx] for xx in feats_box.active])) &
            (DF["sortingMethod"].isin([SORTS[xx] for xx in sort_box.active])) 
           ].copy()
    

    source = ColumnDataSource(df)         
    
    x_title = x.value.title()
    y_title = y.value.title()
    
    kw = dict()
    if x.value in discrete:
        kw['x_range'] = sorted(set(xs))
    if y.value in discrete:
        kw['y_range'] = sorted(set(ys))
    kw['title'] = "%s vs %s" % (x_title, y_title)
    
    
    
    
    if x.value in discrete:
        p.xaxis.major_label_orientation = pd.np.pi / 4
    
    #sz = 9
    if size.value != 'None' and not size.value in discrete :
        if len(set(df[size.value])) > N_SIZES:
            groups = pd.qcut(df[size.value].values, N_SIZES, duplicates='drop')
        else:
            groups = pd.Categorical(df[size.value])
        sz = [SIZES[xx] for xx in groups.codes]
    elif size.value != 'None' :
        groups = pd.Categorical(df[size.value])
        sz = [SIZES[xx] for xx in groups.codes]
    
    n_size=len(np.unique(groups.codes))
    s_size=N_SIZES//n_size
    df["sz"]=[SIZES[xx*s_size] for xx in df[size.value].astype("category").cat.codes]
    source = ColumnDataSource(df)
    #c = "#31AADE"
    if color.value != 'None' and not color.value in discrete :
        if len(set(df[color.value])) > N_COLORS:
            groups = pd.qcut(df[color.value].values, N_COLORS, duplicates='drop')
        else:
            groups = pd.Categorical(df[color.value])
        c = [COLORS[xx] for xx in groups.codes]
    elif color.value != 'None' :
        groups = pd.Categorical(df[color.value])
    
    n_col=len(np.unique(groups.codes))
    s_col=N_COLORS//n_col
    c=factor_cmap(color.value, [COLORS[i*s_col] for i in range(n_col)], DF[color.value].unique())
    
    hover = HoverTool(tooltips=[
    ("RNA type", "@rna_type"),
    ("Features", "@features"),
    ("Grouping Method", "@group"),
    ("Sorting Method", "@sortingMethod"),
    (x.value, "".join(["@",x.value])),
    (y.value, "".join(["@",y.value])),
    ])
    p = figure(height=600, width=800, tools=['pan,box_zoom,reset,wheel_zoom']+[hover], **kw)
    p.xaxis.axis_label = x_title
    p.yaxis.axis_label = y_title
    p.circle(x=x.value,
             y=y.value,
             color=c,
             size="sz",
             source=source,
             line_color="white",
             alpha=0.6,
             hover_color='white',
             hover_alpha=0.5,
             legend_field=color.value
    )
    p.legend.location = 'bottom_right'
    return p


def update(attr, old, new):
    layout.children[1] = create_figure()


x = Select(title='X-Axis', value='Delta_mean', options=variable)
x.on_change('value', update)

y = Select(title='Y-Axis', value='all_max', options=variable)
y.on_change('value', update)

size = Select(title='Size', value='rna_type', options=[None] + discrete)
size.on_change('value', update)

color = Select(title='Color', value='sortingMethod', options=discrete)
color.on_change('value', update)

marker = Select(title='Marker', value=None, options=[None]+discrete)
marker.on_change('value', update)

div_rna = Div(text="""<b>RNA type</b>""")
rna_box = CheckboxGroup(labels = RNAS, active = list(range(len(RNAS))))
rna_box.on_change('active', update)

div_group = Div(text="""<b>Genes groupings</b>""")
group_box = CheckboxGroup(labels = GROUPS, active = list(range(len(GROUPS))))
group_box.on_change('active', update)

div_feats = Div(text="""<b>Features groups</b>""")
feats_box = CheckboxGroup(labels = FEATS, active = list(range(len(FEATS))))
feats_box.on_change('active', update)

div_sort = Div(text="""<b>Sorting Method</b>""")
sort_box = CheckboxGroup(labels = SORTS, active = list(range(len(SORTS))))
sort_box.on_change('active', update)


controls = column(x, y, color, size,
                  div_rna, rna_box,
                  div_group, group_box,
                  div_feats, feats_box,
                  div_sort, sort_box,
                  width=200)
layout = row(controls, create_figure())

curdoc().add_root(layout)
curdoc().title = "CNAM CFDNA"
#output_file(filename="custom_filename.html", title="Static HTML file")
