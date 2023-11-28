import pandas as pd
import matplotlib.pyplot as plt


def convertxvg2png(xvg_file):
    def check_if_value_found(value):
        if value:
            return value[0]
        else:
            return ''


    with open(xvg_file) as inp:
        data = inp.readlines()

    data_coord = [i for i in data if not (i.startswith('#') or i.startswith('@'))]
    title = check_if_value_found([i for i in data if ('@' in i and 'title' in i)]).replace('@    title ','').replace('"','')
    subtitle = check_if_value_found([i for i in data if ('@' in i and 'subtitle' in i)]).replace('@ subtitle ','').replace('"','')
    legend = [i for i in data if ('@ s' in i and 'legend' in i)]

    xaxis = [i for i in data if (i.startswith('@') and 'xaxis' in i)]
    yaxis = [i for i in data if (i.startswith('@') and 'yaxis' in i)]

    coords = [[float(j) for j in i.strip().split(' ') if j] for i in data_coord]

    xaxis = xaxis[0].strip().replace('@    xaxis  label ', '').replace('"', '') if xaxis else 'OX'
    yaxis = yaxis[0].strip().replace('@    yaxis  label ', '').replace('"', '') if yaxis else 'OY'

    plot1 = None
    plt.ioff()
    plt.rcParams.update({'font.size': 15})
    plt.figure(figsize=(15, 12))
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)

    if coords and len(coords[0]) == 2:
        d = pd.DataFrame(coords, columns=[xaxis, yaxis])
        plot1 = plt.plot(d[xaxis], d[yaxis], marker='o', linewidth=2, markersize=5)

    elif coords and len(coords[0]) > 2 and legend:
        d = pd.DataFrame(coords, columns=[xaxis]+legend)
        plt.legend(legend, loc='upper right', bbox_to_anchor=(0, 0), borderaxespad=-1)
        plot1 = plt.plot(d[xaxis], d[legend], marker='o', linewidth=2, markersize=5)
    if plot1:
        plot1 = plot1[0]
        plt.legend(legend)
        plot1 = plot1.figure.suptitle(subtitle)
        plot1.figure.savefig(xvg_file.replace('.xvg','.png'), bbox_inches="tight")