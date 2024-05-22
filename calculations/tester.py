import os.path
import pandas as pd
import matplotlib.pyplot as plt
from numpy import linspace, radians, sin, cos, mean
from tools import to_roman, isnum, export2file, rounding
from colorama import Fore


def Figures(type_fig, *args, **kwargs) -> tuple[list[float]]:
    x0 = kwargs.get('x0', 0)
    y0 = kwargs.get('y0', 0)
    x, y = [], []

    if type_fig.strip().lower() in ('inlet'):
        x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4]
        y = [y0 + 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4]
    elif type_fig.strip().lower() in ('compressor'):
        x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4, x0 - 0.4]
        y = [y0 + 0.4, y0 + 0.2, y0 - 0.2, y0 - 0.4, y0 + 0.4]
    elif type_fig.strip().lower() in ('combustionchamber'):
        x = [0.4 * cos(alpha) + x0 for alpha in linspace(0, radians(360), 360)]
        y = [0.4 * sin(alpha) + y0 for alpha in linspace(0, radians(360), 360)]
    elif type_fig.strip().lower() in ('turbine'):
        x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4, x0 - 0.4]
        y = [y0 + 0.2, y0 + 0.4, y0 - 0.4, y0 - 0.2, y0 + 0.2]
    elif type_fig.strip().lower() in ('nozzle',):
        x = [x0 - 0.4, x0, x0 + 0.4, x0 + 0.4, x0, x0 - 0.4, x0 - 0.4]
        y = [y0 + 0.4, y0 + 0.2, y0 + 0.4, y0 - 0.4, y0 - 0.2, y0 - 0.4, y0 + 0.4]
    elif type_fig.strip().lower() in ('outlet',):
        x = [x0 + 0.4, x0 - 0.4, x0 - 0.4, x0 + 0.4]
        y = [y0 + 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4]
    elif type_fig.strip().lower() in ('heatexchanger',):
        x = [x0 - 0.4, x0 + 0.4, x0 + 0.4, x0 - 0.4, x0 - 0.4]
        y = [y0 + 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4, y0 + 0.4]
    elif type_fig.strip().lower() in ('load',):
        x = [x0 - 0.4, x0, x0 + 0.4, x0 - 0.4]
        y = [y0 - 0.4, y0 + 0.4, y0 - 0.4, y0 - 0.4]
    return x, y


def input_keys(df) -> dict:
    while True:
        z = input('z = ').strip()
        if z in df:
            break
        else:
            print(f'{Fore.RED}There is no such parameter!')
            print(f'{Fore.YELLOW}Possible options:', df.columns.to_list())

    while True:
        x = input('x = ').strip()
        if x in df or not x:
            break
        else:
            print(f'{Fore.RED}There is no such parameter!')
            print(f'{Fore.YELLOW}Possible options:', f'"" {Fore.YELLOW}or', df.columns.to_list())

    while True:
        y = input('y = ').strip()
        if y in df or not y:
            break
        else:
            print(f'{Fore.RED}There is no such parameter!')
            print(f'{Fore.YELLOW}Possible options:', f'"" {Fore.YELLOW}or', df.columns.to_list())

    if not x and not y: return {'z': z}
    if not y: return {'z': z, 'x': x}
    return {'z': z, 'x': x, 'y': y}


def show(df, keys, savefig=False):
    res = dict()
    suptitle = 'gte_analysis'
    plt.figure(figsize=(9, 9))
    plt.suptitle(suptitle, fontsize=16, fontweight='bold')

    D = len(keys)  # мерность пространства
    if 1 == D:
        d = df.to_dict('list')
        # hist
        title = f"{keys['z']}"
        plt.title(title, fontsize=14)
        plt.grid(True)
        plt.hist(d[keys['z']], label=f"{rounding(mean(d[keys['z']]), 4)}")
        plt.xlabel(keys['z'], fontsize=14)
        plt.ylabel('density []', fontsize=14)
    elif 2 == D:
        d = df.sort_values(by=[keys['x']], ascending=True).to_dict('list')
        for i in range(len(df.index)):
            key = ''
            for k in d:
                if k not in keys.values(): key += k + ' = ' + str(rounding(d[k][i], 2)) + '; '
            if key in res:
                res[key][0].append(d[keys['x']][i])
                res[key][1].append(d[keys['z']][i])
            else:
                res[key] = [[d[keys['x']][i]], [d[keys['z']][i]]]

        title = f"{keys['z']}({keys['x']})"
        plt.title(title, fontsize=14)
        plt.grid(True)
        for label, xyz in res.items(): plt.plot(*xyz,
                                                linestyle='-', linewidth=1.5,
                                                marker='o', markersize=3,
                                                label=label)
        plt.xlabel(keys['x'], fontsize=14)
        plt.ylabel(keys['z'], fontsize=14)
    elif 3 == D:
        d = df.sort_values(by=[keys['x'], keys['y']], ascending=True).to_dict('list')
        for i in range(len(df.index)):
            key = ''
            for k in d:
                if k not in keys.values(): key += k + ' = ' + str(rounding(d[k][i], 2)) + '; '
            if key in res:
                res[key][0].append(d[keys['x']][i])
                res[key][1].append(d[keys['y']][i])
                res[key][2].append(d[keys['z']][i])
            else:
                res[key] = [[d[keys['x']][i]], [d[keys['y']][i]], [d[keys['z']][i]]]

        title = f"{keys['z']}({keys['x']}, {keys['y']})"
        plt.title(title, fontsize=14)
        for label, xyz in res.items(): plt.contourf(*xyz, levels=30, cmap='plasma',
                                                    label=label)  # cmap='vinidis'
        plt.colobar(label=keys['z'])

    else:
        raise '4D dimention!'

    for k, v in res.items(): print(k, v)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    if savefig: export2file(plt, file_name=suptitle, file_type='png', show_time=True)
    plt.show()


def analyse(df):
    while True:
        show(df, input_keys(df), savefig=True)


if __name__ == '__main__':
    pass



