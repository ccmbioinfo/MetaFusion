#!/usr/bin/env python
import pandas as pd
import plotly.graph_objs as go
#from plotly.offline import iplot, init_notebook_mode
#import plotly.io as pio
import plotly as plot
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('category_count_file', action='store', help = 'a file containing the counts for each gene category, per caller')
parser.add_argument('output_dir', action = 'store', help = 'the output directory within which the plot file will be generated')
args = parser.parse_args()

# define input file
category_count_file = args.category_count_file

# load data frame
df1 = pd.read_csv(category_count_file, sep='\t', names=None)
df1 = df1.set_index('tool')
df1 = df1.transpose()

layout = dict(title = 'Category counts for 5 fusion callers',
              yaxis = dict(title='Number of Called Fusions'),
              xaxis = dict(title='Gene Fusion Category')
             )

data = [go.Bar(x=df1[tool].index, y=df1[tool], name=tool) for tool in list(df1.columns.values) ]
fig = dict(data=data, layout=layout)
plot.offline.plot(fig, filename=os.path.join(args.output_dir, 'barplot_fusion_categories.html'))
