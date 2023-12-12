import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import os
from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import flux_variability_analysis
import numpy as np


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
    # sns.scatterplot(x='tsne 1', y='tsne 2', hue='factor', data=data, ax=ax, s=60)

    sns.scatterplot(x='tsne 1', y='tsne 2', hue='factor', data=data, ax=ax, s=60,
                    palette=dict(green="#2ca02c",
                                 mature="#d62728"))

    # a = pd.concat({'x': data['tsne 1'], 'y': data['tsne 2']}, axis=1)
    # for i, point in a.iterrows():
    #     ax.text(point['x'] + .02, point['y'], str(i), fontdict={'size': 5})

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

    simulations_folder = 'C:/Users/BiSBII/Documents/MM_ML/data/model_simulations'

    all_models = os.listdir(models_folder)

    for modelfile in all_models:
        specific_model = read_sbml_model(os.path.join(models_folder, modelfile))
        res = specific_model.optimize()

        reactions = [r for r in res.fluxes.index if not r.startswith('EX_')]
        fluxes = res.fluxes[reactions]

        name = '_'.join(modelfile.split('_')[:3])

        fluxes.name = name

        fluxes.to_csv(os.path.join(simulations_folder, name + '.csv'))


def split_reversible_reactions(model_to_sample):
    # exchanges_demands_sinks = [reaction.id for reaction in model_to_sample.exchanges] + [reaction.id for reaction in model_to_sample.demands] + [reaction.id for reaction in model_to_sample.sinks]
    # exchanges_demands_sinks = set(exchanges_demands_sinks)
    drains = [r for r in model_to_sample.reactions if r.id.startswith('EX_')]
    new_reactions = []
    for reaction in model_to_sample.reactions:
        if reaction not in drains:
            if reaction.lower_bound < 0 < reaction.upper_bound:
                new_reaction = reaction.copy()
                new_reaction.id = reaction.id + "_reverse"
                new_reaction.lower_bound = 0
                new_reaction.upper_bound = -reaction.lower_bound
                for metabolite, coefficient in new_reaction.metabolites.items():
                    new_reaction.add_metabolites({metabolite: -coefficient})
                    new_reaction.add_metabolites({metabolite: -coefficient})
                new_reactions.append(new_reaction)
                reaction.lower_bound = 0
    model_to_sample.add_reactions(new_reactions)
    return model_to_sample


def model_fva():
    models_folder = ('C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2023/results_troppo/ALL_BERRY/'
                     'reconstructed_models')

    simulations_folder = 'C:/Users/BiSBII/Documents/MM_ML/data/model_fc'

    all_models = os.listdir(models_folder)

    for modelfile in all_models[5:]:
        specific_model = read_sbml_model(os.path.join(models_folder, modelfile))

        specific_model_irr = split_reversible_reactions(specific_model)

        res = flux_variability_analysis(specific_model_irr, fraction_of_optimum=0.8)

        reactions = [r for r in res.index if not r.startswith('EX_')]
        flux_max = res.loc[reactions, 'maximum']
        flux_min = res.loc[reactions, 'minimum']

        fc = flux_max.subtract(flux_min)

        name = '_'.join(modelfile.split('_')[:3])

        fc.name = name

        fc.to_csv(os.path.join(simulations_folder, name + '.csv'))


def create_fluxomics_df():
    simulations_folder = 'C:/Users/BiSBII/Documents/MM_ML/data/model_fc'

    df_files = os.listdir(simulations_folder)

    all_dfs = []
    for f in df_files:
        df = pd.read_csv(os.path.join(simulations_folder, f), index_col=0)
        all_dfs.append(df)

    df_fluxomics = pd.concat(all_dfs, axis=1)
    df_fluxomics.to_csv('C:/Users/BiSBII/Documents/MM_ML/data/fluxomics_fc.csv')


def get_correlations():
    folder = 'C:/Users/BiSBII/Documents/MM_ML/data/DIABLO_INPUT'

    df = pd.read_excel(os.path.join(folder, 'correlations_python.xlsx'), index_col=0)

    np.fill_diagonal(df.values, 0)

    to_keep = [c for c in df.columns if not c.startswith('Vit') and '__' in c]

    genes = [i for i in df.index if i.startswith('Vit')]

    df_mets = df.loc[genes, to_keep]
    df_mets_08 = df_mets[df_mets.columns[(abs(df_mets) > 0.95).any()]]

    for col in df_mets_08:

        neg = df_mets_08.loc[df_mets_08[col] < -0.95, col].sort_values()

        pos = df_mets_08.loc[df_mets_08[col] > 0.95, col].sort_values(ascending=False)

        colname = col.replace('.', '_')
        with pd.ExcelWriter(os.path.join(folder, 'correlations_all_python', colname + '.xlsx')) as writer:
            neg.to_excel(writer, sheet_name='negative_correlation')
            pos.to_excel(writer, sheet_name='postive_correlation')


def get_paths():
    model_file = 'C:/Users/BiSBII/Documents/plantdb/reconstruction_results/vvinif2023/vvinif2023_FINAL.xml'
    model = read_sbml_model(model_file)
    model_genes = {}
    for gene in model.genes:
        new_id = gene.id.split('_')[0]
        if new_id not in model_genes:
            model_genes[new_id] = [gene]
        else:
            model_genes[new_id].append(gene)

    # folder = 'C:/Users/BiSBII/Documents/MM_ML/data/DIABLO_INPUT/correlations'
    folder = 'C:/Users/BiSBII/Documents/MM_ML/data/DIABLO_INPUT'
    df = pd.read_csv(os.path.join(folder, 'genes_loadings.csv'), header=None)
    genes = df[0].to_list()
    # files = [f for f in os.listdir(folder) if f.endswith('.xlsx')]

    # for met_file in files:
    #     neg = pd.read_excel(os.path.join(folder, met_file), index_col=0, sheet_name='negative_correlation')
    #     pos = pd.read_excel(os.path.join(folder, met_file), index_col=0, sheet_name='postive_correlation')

    neg_paths = []
    for g in genes:
    # for g in neg.index:
        if g.startswith('Vit'):
            objc = model_genes[g][0]
            all_paths = []
            for reac in objc.reactions:
                groups = [(c.id, c.name) for c in model.get_associated_groups(reac)]
                for group in groups:
                    if group not in all_paths:
                        all_paths.append(group)
            neg_paths.append(all_paths)
        else:
            neg_paths.append([])

    for i in range(len(genes)):
        print(genes[i], ':', neg_paths[i])

        # pos_paths = []
        # for g in pos.index:
        #     if g.startswith('Vit'):
        #         objc = model_genes[g][0]
        #         all_paths = []
        #         for reac in objc.reactions:
        #             groups = [(c.id, c.name) for c in model.get_associated_groups(reac)]
        #             for group in groups:
        #                 if group not in all_paths:
        #                     all_paths.append(group)
        #         pos_paths.append(all_paths)
        #     else:
        #         pos_paths.append([])

        # neg['paths'] = neg_paths
        # pos['paths'] = pos_paths
        #
        # with pd.ExcelWriter(os.path.join(folder, 'paths', met_file.replace('.xlsx', '') + '_path.xlsx')) as writer:
        #     neg.to_excel(writer, sheet_name='negative_correlation')
        #     pos.to_excel(writer, sheet_name='positive_correlation')



if __name__ == '__main__':
    # add_drains_transporters_models()
    # model_simulation()
    # model_fva()
    # create_fluxomics_df()
    get_correlations()
    # get_paths()
