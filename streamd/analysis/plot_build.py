import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib import rcParams
import plotly.offline
import plotly.express as px

# figure size in inches
#rcParams['figure.figsize'] = 15,10

#sns.set_theme(rc={'figure.figsize':(11.7,8.27)})
sns.set_context("paper", rc={"font.size":15,"axes.titlesize":15,"axes.labelsize":15},
                font_scale=1.5)

plt.ioff()

def plot_rmsd(rmsd_df, system_name, out):
    plot = rmsd_df.set_index('time(ns)').plot(title=f"RMSD of {system_name}")
    plt.ylabel("RMSD (Å)")
    plt.xlabel("Time (ns)")
    plt.legend(loc='lower right', ncol=len(rmsd_df.columns) // 2, frameon=False)
    plot.figure.savefig(out, dpi=350, bbox_inches="tight")
    plt.clf()
    plt.close('all')

def plot_rmsd_mean_std(data, paint_by_col, show_legend, out_name, title=None):
    #pd.DataFrame.iteritems = pd.DataFrame.items
    # df = pd.read_csv(rmsd_mean_std_fname, sep='\t')
    # g = sns.FacetGrid(df, row='rmsd_system', col='time',
    #                   height=8, aspect=2,
    #                   hue='rmsd_system',
    #                   margin_titles=True)
    # g.map(sns.scatterplot, 'RMSD_mean', 'RMSD_std')
    # g.set_axis_labels('RMSD Mean', 'RMSD Std')
    # g.set_titles(col_template="{col_name}", row_template="{row_name}")
    # g.refline(x=2, y=0.5, color='grey')
    # plt.annotate(f"optimal zone", xy=(2, 0.5),
    #              xycoords=plt.gca().get_yaxis_transform(), ha="right")
    # g.tight_layout()
    # g.savefig(out+'.png', dpi=350)
    hover_data = {'system': True, 'time_range': False, 'rmsd_system': False}
    if 'ligand_name' in data:
        hover_data['ligand_name'] = True
    if 'directory' in data:
        hover_data['directory'] = True

    fig = px.scatter(
        data,
        x='RMSD_mean',
        y='RMSD_std',
        color=paint_by_col,
        facet_row='rmsd_system',
        facet_col='time_range',
        title= title if title else 'RMSD Mean vs RMSD Std',
        labels={'RMSD_mean': 'RMSD Mean', 'RMSD_std': 'RMSD Std'},
        hover_data=hover_data,
    )

    fig.add_hline(y=0.5, line_width=0.5, line_dash="dash", row='all', col='all')
    fig.add_vline(x=5, line_width=0.5, line_dash="dash", row='all', col='all')

    fig.add_shape(type="rect",
                  x0=0, y0=0, x1=5, y1=0.5,
                  line=dict(
                      color="darkgreen",
                      width=2,
                      dash='dash'
                  ),
                  fillcolor='rgba(0, 255, 0, 0.2)',
                  opacity=0.4,
                  row='all', col='all'
                  )

    fig.update_annotations(font_size=18)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.update_layout(
        margin=dict(l=0, r=0, b=0, t=80),
        plot_bgcolor='white',
        height=650,
        width=1400,
        title_x=0.5,
        showlegend=show_legend
     )
    fig.update_xaxes(showline=True, linecolor='black', linewidth=0.1)
                     #rangemode='nonnegative')
    fig.update_yaxes(showline=True, linecolor='black', linewidth=0.1, mirror=True)
    fig.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True))
    fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
    plotly.offline.plot(fig, filename=out_name,  auto_open=False)