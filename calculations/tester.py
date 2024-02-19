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


def correlation(data, *args, method='pearson', only_nums=False, **kwargs):
    if method.strip().lower() in ('pearson', 'kendall', 'spearman'):
        cor = data[1:].corr(method=method)
    else:
        print('method must be in "pearson" or "kendall" or "spearman"!')
        return
    if only_nums:  # удаление всех строк и столбцов, содержащие только nan
        cor.dropna(how='all', axis=0, inplace=True)  # удаление всех строк содержащих только nan
        cor.dropna(how='all', axis=1, inplace=True)  # удаление всех столбцов содержащих только nan
    return cor


def show_correlation(df, show_num=True, cmap='bwr', units='', rnd=4, savefig=False, **kwargs):
    cor = correlation(df, **kwargs)
    if cor is None: return
    if not isnum(rnd, type_num='int') or rnd < 0: rnd = 4
    if units == '%': cor = cor * 100
    cor = cor.round(rnd)

    fig, ax = plt.subplots(figsize=(len(df.columns) / 2.54, len(df.columns) / 2.54))
    ax.set_title('Correlation matrix | Матрица корреляции', fontsize=16, fontweight='bold')
    ax.set_aspect('equal')
    im = ax.imshow(cor, interpolation='nearest',
                   cmap=cmap)  # RGB: 'jet', 'turbo', # blue vs red: 'bwr', # 2side: 'twilight','twilight_shifted'

    cbar = fig.colorbar(im, orientation='vertical')
    if units == '':
        cbar.set_label('Color Intensity []', fontsize=14)
        cbar.set_ticks(linspace(-1, 1, 21))
    else:
        cbar.set_label('Color Intensity [%]', fontsize=14)
        cbar.set_ticks(linspace(-100, 100, 21))

    ax.set_xticks(range(len(cor.columns)), cor.columns, fontsize=14, rotation=90)
    ax.set_yticks(range(len(cor.columns)), cor.columns, fontsize=14, rotation=0)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    if show_num:
        for row in range(len(cor.columns)):
            for col in range(len(cor.columns)):
                ax.text(row, col, cor.to_numpy()[row, col],
                        ha="center", va="center", color='black', fontsize=12, rotation=45)

    if savefig: export2file(plt, file_name='correlation', show_time=True, file_type='png', **kwargs)

    plt.show()


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
    file = 'exports/GTE АИ-222-25.pkl'

    if os.path.isfile(file):
        file_name, file_extension = os.path.splitext(file)
        if file_extension not in ('.pkl', '.csv', '.xlsx'):
            print('Неизвестное расширение!')
    else:
        print('Такого файла не существет!')
        exit()

    if file_extension == '.pkl':
        raw_df = pd.read_pickle(file)
    elif file_extension == '.csv':
        raw_df = pd.read_csv(file)
    elif file_extension == '.xlsx':
        raw_df = pd.read_excel(file)

    print(raw_df)
    print(raw_df.info(memory_usage='deep'))

    df = raw_df
    data = raw_df.to_dict('list')
    print(data)

    scheme = dict()
    c = 0
    for k in data:
        if 'contour' in k:
            c += 1
            scheme[c] = data[k][0].split('+')
    del c

    fg = plt.figure(figsize=(13, 6))  # размер в дюймах
    fg.suptitle('GTE scheme', fontsize=14, fontweight='bold')
    gs = fg.add_gridspec(1, len(scheme))  # строки, столбцы

    for contour in scheme:
        fg.add_subplot(gs[0, contour - 1])
        plt.grid(True)
        plt.axis('square')
        plt.title('contour ' + to_roman(contour) + ' | ' + 'контур ' + to_roman(contour), fontsize=14)
        plt.xlim(0, len(scheme[contour]))
        plt.ylim(0, 1)
        plt.xticks(linspace(0, len(scheme[contour]), len(scheme[contour]) + 1))
        plt.yticks(linspace(0, 1, 1 + 1))

        x0, y0 = 0.5, 0.5

        for node in scheme[contour]:
            # print(node, *Figures(node, x0=x0, y0=y0))
            plt.plot(*Figures(node, x0=x0, y0=y0), color='black', linewidth=3)
            x0 += 1

    plt.show()

    list2del = []
    for k in data:
        for v in data[k]:
            if not isnum(v):
                list2del.append(k)
                break
    for i, k in enumerate(list2del): del data[k]
    del list2del
    print(data)

    clean_df = pd.DataFrame(data)
    print(clean_df)

    cor = correlation(clean_df, method='pearson', only_nums=True)
    print(cor)

    show_correlation(clean_df, only_nums=False, show_num=False, units='%', savefig=False)
    show_correlation(clean_df, only_nums=True, show_num=True, units='', rnd=2, savefig=True, dpi=300)

    analyse(df)
