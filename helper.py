import numpy as np
import pandas as pd
from sklearn.metrics import zero_one_loss
from sklearn.model_selection import StratifiedKFold
from scipy import stats
from scipy.cluster import hierarchy
from statsmodels.stats.multicomp import pairwise_tukeyhsd

ALPHA = .05


def plot_dendrogram(model, **kwargs):

    # Children of hierarchical clustering
    children = model.children_

    # Distances between each pair of children
    # Since we don't have this information, we can use a uniform one for plotting
    distance = np.arange(children.shape[0])

    # The number of observations contained in each cluster level
    no_of_observations = np.arange(2, children.shape[0]+2)

    # Create linkage matrix and then plot the dendrogram
    linkage_matrix = np.column_stack([children, distance, no_of_observations]).astype(float)

    # Plot the corresponding dendrogram
    hierarchy.dendrogram(linkage_matrix, **kwargs)


def evaluate_classifier(classifier, X, y):
    k_fold = StratifiedKFold(n_splits=3)

    losses = []
    for train_index, test_index in k_fold.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        loss = zero_one_loss(y_test, classifier.fit(X_train, y_train).predict(X_test))
        print(f'\tDataset sizes: ({train_index.shape[0]}, {test_index.shape[0]}), loss: {loss:.3f}')
        losses.append(loss)
    print(f'k-Fold mean loss: {np.mean(losses):.3f}, std: {np.std(losses):.3f}')

    classifier.fit(X, y)
    loss = zero_one_loss(y, classifier.predict(X))
    print(f'Loss of training on full dataset: {loss:.2f}')

    print('-'*15, '\n')
    return classifier


def split_into_groups(data, gene_name):
    group_names, groups = [], []
    for key, group in data.groupby('group'):
        groups.append(group[gene_name].values)
        group_names.append(key)
    return group_names, groups


def apply_anova(data):
    p_values = []
    genes = data.drop(columns=['group']).columns
    for i, col in enumerate(genes):
        group_names, groups = split_into_groups(data, col)
        res = stats.f_oneway(*groups)

        shapiro_res = stats.shapiro(np.concatenate(groups))  # normality test
        levene_res = stats.levene(*groups, center='mean')  # homodestacity test

        p_values.append((col, res.pvalue, shapiro_res[1], levene_res.pvalue))

        if i % 100 == 0:
            print('Progress {:2.0%}'.format((i / (genes.shape[0]))), end='\r')

    anova_table = pd.DataFrame(p_values, columns=['gene', 'p_value', 'shapiro_p_value',
                                                  'levene_p_value'])
    print('Found {} genes that influence'.format((anova_table['p_value'] < ALPHA).sum()),
          'health conditions according to ANOVA tests.')
    print('Found {} genes that influence'.format(((anova_table['p_value'] < ALPHA) &
                                                 (anova_table['shapiro_p_value'] > ALPHA) &
                                                 (anova_table['levene_p_value'] > ALPHA)).sum()),
          'health conditions according to ANOVA Shapiro-Wilks and Levene tests.')

    return anova_table


def apply_tukey_hsd(data, anova_table):
    p_values = []
    genes = anova_table.loc[(anova_table['p_value'] < ALPHA) &
                            (anova_table['shapiro_p_value'] > ALPHA) &
                            (anova_table['levene_p_value'] > ALPHA),
                            'gene']
    for i, col in enumerate(genes):
        res = pairwise_tukeyhsd(data[col], data['group'])
        p_values.append((col,) + tuple(res.pvalues))

        if i % 100 == 0:
            print('Progress {:2.0%}'.format((i / (genes.shape[0]))), end='\r')

    col_names = [f'{res.groupsunique[i]}-{res.groupsunique[j]}'
                 for i in range(res.groupsunique.shape[0])
                 for j in range(i+1, res.groupsunique.shape[0])]
    tukey_table = pd.DataFrame(p_values, columns=['gene'] + col_names)

    return tukey_table
