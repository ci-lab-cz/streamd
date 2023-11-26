import seaborn
import pandas as pd
import matplotlib.pyplot as plt
plt.ioff()

def convertxvg2png(xvg_file):
    def check_if_value_found(value):
        if value:
            return value[0]
        else:
            return ''


    with open(xvg_file) as inp:
        data = inp.readlines()

    data_coord = [i for i in data if not (i.startswith('#') or i.startswith('@'))]
    title = check_if_value_found([i for i in data if ('@' in i and 'title' in i)])
    subtitle = check_if_value_found([i for i in data if ('@' in i and 'subtitle' in i)])
    legend = [i for i in data if ('@ s' in i and 'legend' in i)]

    xaxis = [i for i in data if (i.startswith('@') and 'xaxis' in i)]
    yaxis = [i for i in data if (i.startswith('@') and 'yaxis' in i)]

    coords = [[float(j) for j in i.strip().split(' ') if j] for i in data_coord]

    xaxis = xaxis[0].strip().replace('@    xaxis  label ', '').replace('"', '') if xaxis else 'OX'
    yaxis = yaxis[0].strip().replace('@    yaxis  label ', '').replace('"', '') if yaxis else 'OY'

    plot = None
    plt.figure(figsize=(15, 8))
    plt.legend(loc='upper right', bbox_to_anchor=(0, 0), borderaxespad=-1)

    if coords and len(coords[0]) == 2:
        d = pd.DataFrame(coords, columns=[xaxis, yaxis])
        plot = seaborn.lineplot(d, x=xaxis, y=yaxis, marker="o").set_title(title)

    elif coords and len(coords[0]) > 2 and legend:
        d = pd.DataFrame(coords, columns=[xaxis]+legend).pivot_table(index=xaxis, values=legend)
        plot = seaborn.lineplot(d, marker="o").set_title(title)

    if plot:
        plot = plot.figure.suptitle(subtitle)
        plot.figure.savefig(xvg_file.replace('.xvg','.png'), bbox_inches="tight")

