import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import os
from cobra.io import read_sbml_model, write_sbml_model

def run_tsne(n_components: int, data: pd.DataFrame):
    tsne = TSNE(n_components, random_state=42, perplexity=30)
    tsne_result = tsne.fit_transform(data)

    columns = [f'tsne {i + 1}' for i in range(n_components)]

    df_tsne = pd.DataFrame(data=tsne_result, index=data.index, columns=columns)

    return df_tsne


def plot_tsne(data: pd.DataFrame, name_fig: str, title: str):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    sns.set_style("darkgrid")
    sns.scatterplot(x='tsne 1', y='tsne 2', hue='factor', data=data, ax=ax, s=60)

    # a = pd.concat({'x': data['tsne 1'], 'y': data['tsne 2']}, axis=1)
    # for i, point in a.iterrows():
    #     ax.text(point['x'] + .01, point['y'], str(i), fontsize='x-small')

    ax.set_title(title, fontsize=20)
    ax.legend(loc='lower left', borderaxespad=0.0)
    fig_path = os.path.join('C:/Users/BiSBII/Documents/MM_ML', name_fig + '.png')
    plt.savefig(fig_path)
    plt.show()


def add_drains_transporters_models():
    models_folder = ('C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2023/results_troppo/ALL_BERRY/'
                     'reconstructed_models')

    final_model = 'C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2023/vvinif2023_FINAL.xml'

    model = read_sbml_model(final_model)

    drains = [r for r in model.reactions if r.id.startswith('EX_')]

    transporters_t = [r for r in model.reactions if r.id.startswith('T_')]

    all_models = os.listdir(models_folder)

    for modelfile in all_models:

        specific_model = read_sbml_model(os.path.join(models_folder, modelfile))

        drains_to_add = []
        transporters_to_add = []

        for d in drains:
            try:
                specific_model.reactions.get_by_id(d.id)
            except KeyError:
                drains_to_add.append(d)

        for t in transporters_t:
            try:
                specific_model.reactions.get_by_id(t.id)
            except KeyError:
                transporters_to_add.append(t)

        specific_model.add_reactions(drains_to_add + transporters_to_add)

        g_drains = specific_model.groups.get_by_id('drains')
        g_drains.add_members(drains_to_add)

        g_transp = specific_model.groups.get_by_id('transporters')
        g_transp.add_members(transporters_to_add)

        write_sbml_model(specific_model, os.path.join(models_folder, modelfile))


def model_simulation():
    models_folder = ('C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2023/results_troppo/ALL_BERRY/'
                     'reconstructed_models')

    all_models = os.listdir(models_folder)

    for modelfile in all_models:
        specific_model = read_sbml_model(os.path.join(models_folder, modelfile))
        res = specific_model.optimize()
        fluxes = res.fluxes
        reactions = res.fluxes.index
        print(reactions)

        break


if __name__ == '__main__':
    add_drains_transporters_models()
    model_simulation()