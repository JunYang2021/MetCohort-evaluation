# fig 6, figs17-19

import re

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os
import scipy.stats as stats

font = {'family': 'Arial', 'size': 7}
plt.rc('font', **font)

# Fig 6a, Fig s17
"""# Quantitation evaluation in BM21
# 1. Get file names and concentration proportion
file_names = []
concentrations = []
con_proportion = {'1024;1': 1024 / 1025,
                  '1;1024': 1 / 1025,
                  '1;16': 1 / 17,
                  '1;1': 1 / 2,
                  '1;256': 1 / 257,
                  '1;4': 1 / 5,
                  '16;1': 16 / 17,
                  '1;64': 1 / 65,
                  '256;1': 256 / 257,
                  '4;1': 4 / 5,
                  '64;1': 64 / 65}
with open('../BM21/Description.txt', 'r') as file:
    for line in file:
        if line.startswith("SUBJECT_SAMPLE_FACTORS"):
            parts = line.split()
            file_names.append(parts[2])
            pattern = r'Veg\(([^)]+)\)'
            con = re.search(pattern, parts[3]).group(1)
            concentrations.append(con_proportion[con])
# Get concentration correlation from xcms
xcms_cor = []
xcms_features = pd.read_excel('../BM21/xcms_bm21.xlsx')
for _, fea in xcms_features.iterrows():
    xcms_data = fea[file_names]
    if len(np.unique(xcms_data)) == 1:
        print('xcms')
        continue
    cor, _ = pearsonr(xcms_data, concentrations)
    xcms_cor.append(abs(cor))
# Get concentration correlation from asari
asari_cor = []
asari_features = pd.read_excel('../BM21/asari-bm21.xlsx')
for _, fea in asari_features.iterrows():
    asari_data = fea[file_names]
    if len(np.unique(asari_data)) == 1:
        print('asari')
        continue
    cor, _ = pearsonr(asari_data, concentrations)
    asari_cor.append(abs(cor))
# Get concentration correlation from MetCohort
MetCohort_cor = []
MetCohort_features = pd.read_excel('../BM21/peakmat-v3-4261-hilic-veg.xlsx')
for _, fea in MetCohort_features.iterrows():
    MetCohort_data = fea[file_names]
    if len(np.unique(MetCohort_data)) == 1:
        print('MetCohort')
        continue
    cor, _ = pearsonr(MetCohort_data, concentrations)
    MetCohort_cor.append(abs(cor))

# Plot the distribution of correlation from three tools
fig, ax = plt.subplots(figsize=(4, 3.5))
# fig, ax = plt.subplots(figsize=(3.5, 2))
x = np.linspace(0, 1, 200)
xcms_density = gaussian_kde(xcms_cor)
xcms_y = xcms_density(x)
# ax.plot(x, xcms_y, color='#e41a1c', label='XCMS')
plt.hist(xcms_cor,bins=40, color='#e41a1c', label='XCMS', alpha=0.7)
asari_density = gaussian_kde(asari_cor)
asari_y = asari_density(x)
# ax.plot(x, asari_y, color='#377eb8', label='asari')
plt.hist(asari_cor, bins=40, color='#377eb8', label='asari', alpha=0.7)
MetCohort_density = gaussian_kde(MetCohort_cor)
MetCohort_y = MetCohort_density(x)
# ax.plot(x, MetCohort_y, color='#4daf4a', label='MetCohort')
plt.hist(MetCohort_cor, bins=40, color='#4daf4a', label='MetCohort', alpha=0.7)
ax.legend()
plt.xlabel('Pearson |r|')
plt.ylabel('Feature number')
# plt.ylabel('Feature density')
plt.tight_layout()
plt.savefig('D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\\figure7\\bm21-corr_renamed.tiff', dpi=300)
# plt.savefig('D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\\figure7\\bm21-corr-density_renamed.tiff', dpi=300)"""


# Plot feature correlation of different feature tables
print('read_MetCohort...')
MetCohort_file = "D:\\Experiments\\peakmat_evaluation\\XG\\peakmat-v3_8600.csv"
MetCohort_df = pd.read_csv(MetCohort_file)
print('read_xcms...')
xcms_file = "D:\\Experiments\\peakmat_evaluation\\XG\\xcms_1769files_20240601.csv"
xcms_df = pd.read_csv(xcms_file)
picture_dir = 'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\\figure7'


def format_formula(formula):
    formatted = ''
    i = 0
    while i < len(formula):
        if formula[i].isdigit():
            formatted += f'$_{{{formula[i]}}}$'
        elif formula[i] == '[' and formula[i + 1] == '1' and formula[i + 2] == '3':
            formatted += '$^{13}$'
            i += 3
        else:
            formatted += formula[i]
        i += 1
    return formatted


def plot_MetCohort_cor(peak1_id, peak2_id, peak1_name, peak2_name):
    # Adjusting for 0-based indexing
    peak1_id -= 1
    peak2_id -= 1

    # Extracting the data for the specified peaks
    peak1_data = MetCohort_df.iloc[peak1_id, 6:]  # Data starts from the 6th column
    peak2_data = MetCohort_df.iloc[peak2_id, 6:]
    print(len(peak1_data))

    r, _ = pearsonr(peak1_data, peak2_data)

    fig, ax = plt.subplots(figsize=(2, 2))

    plt.scatter(peak1_data, peak2_data, s=1, c='#66c2a5', alpha=0.9)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
    # plt.scatter(peak1_data, peak2_data, s=2.5, c='#fc8d62', alpha=0.9)
    plt.xlabel(format_formula(peak1_name))
    plt.ylabel(format_formula(peak2_name))
    plt.figtext(0.4, 0.75, f'MetCohort\nR = {r:.4f}', ha='center', va='center')
    plt.tight_layout()

    picture_name = os.path.join(picture_dir, f'MetCohort_{peak1_id + 1}_{peak2_id + 1}.tiff')
    plt.savefig(picture_name, dpi=300)


def plot_xcms_cor(peak1_id, peak2_id, peak1_name, peak2_name):
    # Adjusting for 0-based indexing
    peak1_id -= 1
    peak2_id -= 1

    # Extracting the data for the specified peaks
    peak1_data = xcms_df.iloc[peak1_id, 10:]  # Data starts from the 10th column
    peak2_data = xcms_df.iloc[peak2_id, 10:]
    print(len(peak1_data))

    r, _ = pearsonr(peak1_data, peak2_data)

    fig, ax = plt.subplots(figsize=(2, 2))

    plt.scatter(peak1_data, peak2_data, s=1, c='#fc8d62', alpha=0.9)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
    plt.xlabel(format_formula(peak1_name))
    plt.ylabel(format_formula(peak2_name))
    # plt.ylim(0,1e6)
    plt.figtext(0.4, 0.75, f'XCMS\nR = {r:.4f}', ha='center', va='center')
    plt.tight_layout()

    picture_name = os.path.join(picture_dir, f'xcms_{peak1_id + 1}_{peak2_id + 1}.tiff')
    plt.savefig(picture_name, dpi=300)


def calculate_distance(row1, row2):
    return np.sqrt((row1['mz'] - row2['mz']) ** 2 * 1e5 + (row1['rt(s)'] - row2['rt']) ** 2)


def print_nan_inf(data):
    for column in data.index:
        value = data[column]
        if pd.isna(value) or np.isinf(value):
            print(f'{column}: {value}')


def sanitize_filename(name):
    """ Replace or remove illegal characters in filenames. """
    name = name.replace(":", "_")  # Replace colons with underscores
    name = name.replace(" ", "_")  # Replace spaces with underscores
    name = re.sub(r"[^\w\s-]", "", name)  # Remove other non-alphanumeric characters
    return name


def compare_integration(feature_table):
    compounds = pd.read_excel(feature_table, sheet_name="metabolites")
    samples_true = pd.read_excel(feature_table, sheet_name="samples")

    delta_mz = 0.01
    delta_rt = 10

    MetCohort_det = []
    xcms_det = []

    for index, row in compounds.iterrows():
        print('compounds:', index)

        # Lists to store features and their distances
        MetCohort_candidates = []
        xcms_candidates = []

        # Find candidates in MetCohort_features
        for _, row_p in MetCohort_df.iterrows():
            if abs(row['mz'] - row_p['mz']) <= delta_mz and abs(row['rt(s)'] - row_p['rt']) <= delta_rt:
                distance = calculate_distance(row, row_p)
                MetCohort_candidates.append((row_p['id'], distance))

        # Select the feature with the smallest distance in MetCohort_features
        if MetCohort_candidates:
            MetCohort_det.append(min(MetCohort_candidates, key=lambda x: x[1])[0])
        else:
            MetCohort_det.append('None')

        # Find candidates in xcms_features
        for _, row_x in xcms_df.iterrows():
            if abs(row['mz'] - row_x['mz']) <= delta_mz and abs(row['rt(s)'] - row_x['rt']) <= delta_rt:
                distance = calculate_distance(row, row_x)
                xcms_candidates.append((row_x['id'], distance))

        # Select the feature with the smallest distance in xcms_features
        if xcms_candidates:
            xcms_det.append(min(xcms_candidates, key=lambda x: x[1])[0])
        else:
            xcms_det.append('None')

    compounds['peakmat_detection'] = MetCohort_det
    compounds['xcms_detection'] = xcms_det
    compounds['correlation_coefficient_peakmat'] = None
    compounds['correlation_coefficient_xcms'] = None
    compounds['kendall_correlation_coefficient_peakmat'] = None
    compounds['kendall_correlation_coefficient_xcms'] = None
    compounds['kendall_correlation_coefficient_peakmat_p'] = None
    compounds['kendall_correlation_coefficient_xcms_p'] = None
    files_name = MetCohort_df.columns[6:]
    files_name = [col for col in files_name if col != 'QC-33']      # This file doesn't have integration data

    for index, compound_row in compounds.iterrows():
        print('compounds cor:', index)
        # if compound_row['compounds'] in {'17-Hydroxyprogesterone', 'linolenyl carnitine'}:
        #     continue

        if compound_row['peakmat_detection'] != 'None':
            MetCohort_data = MetCohort_df.loc[MetCohort_df['id'] == compound_row['peakmat_detection'], files_name]
            reference_data = samples_true.loc[samples_true['compounds'] == compound_row['compounds'], files_name]

            if not MetCohort_data.empty and not reference_data.empty and MetCohort_data.shape == reference_data.shape:
                print_nan_inf(MetCohort_data.iloc[0])
                print_nan_inf(reference_data.iloc[0])
                corr1, _ = pearsonr(MetCohort_data.iloc[0], reference_data.iloc[0])
                compounds.at[index, 'correlation_coefficient_peakmat'] = corr1
                tau, p_value = stats.kendalltau(MetCohort_data.iloc[0], reference_data.iloc[0])
                compounds.at[index, 'kendall_correlation_coefficient_peakmat'] = tau
                compounds.at[index, 'kendall_correlation_coefficient_peakmat_p'] = p_value

        if compound_row['xcms_detection'] != 'None':
            xcms_data = xcms_df.loc[xcms_df['id'] == compound_row['xcms_detection'], files_name]
            reference_data = samples_true.loc[samples_true['compounds'] == compound_row['compounds'], files_name]

            if not xcms_data.empty and not reference_data.empty and xcms_data.shape == reference_data.shape:
                print_nan_inf(xcms_data.iloc[0])
                print_nan_inf(reference_data.iloc[0])
                corr2, _ = pearsonr(xcms_data.iloc[0], reference_data.iloc[0])
                compounds.at[index, 'correlation_coefficient_xcms'] = corr2
                tau, p_value = stats.kendalltau(xcms_data.iloc[0], reference_data.iloc[0])
                compounds.at[index, 'kendall_correlation_coefficient_xcms'] = tau
                compounds.at[index, 'kendall_correlation_coefficient_xcms_p'] = p_value

        if compound_row['peakmat_detection'] != 'None' and compound_row['xcms_detection'] != 'None':
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(3.7, 2))
            ax1.scatter(MetCohort_data, reference_data, s=1, c='#66c2a5', alpha=0.9)
            ax1.set_xlabel('MetCohort quantification')
            ax1.set_ylabel('Trace Finder integration')
            ax1.text(0.3, 0.75, f'R = {corr1:.4f}', ha='center', va='center', transform=ax1.transAxes)
            ax2.scatter(xcms_data, reference_data, s=1, c='#fc8d62', alpha=0.9)
            ax2.set_xlabel('XCMS quantification')
            ax2.text(0.3, 0.75, f'R = {corr2:.4f}', ha='center', va='center', transform=ax2.transAxes)
            plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
            plt.tight_layout()
            sanitized_compound_name = sanitize_filename(compound_row['compounds'])
            _path = 'D:\\Metabolomics2022\\文章\\peakmat文章\\peakmat article 20240513\\pictures\\figure7\\compare_with_manual_renamed\\' + f"{sanitized_compound_name}.png"
            plt.savefig(_path, dpi=300)
            plt.close()

    compounds.to_excel('D:\\Experiments\\peakmat_evaluation\\XG\\quantification_correlation_20240612.xlsx')


def plot_correlation(input_file, output_file):
    features = pd.read_excel(input_file)
    MetCohort_cor = features['correlation_coefficient_peakmat']
    xcms_cor = features['correlation_coefficient_xcms']
    x_MetCohort = np.random.rand(len(MetCohort_cor))
    x_xcms = np.random.rand(len(xcms_cor))

    fig, ax = plt.subplots(figsize=(3.5, 2))
    plt.scatter(x_xcms, xcms_cor, c='#e41a1c', label='XCMS', s=1, alpha=0.8)
    plt.scatter(x_MetCohort, MetCohort_cor, c='#377eb8', label='MetCohort', s=1, alpha=0.8)


    plt.axhline(y=0.8, color='#4daf4a', linestyle='--')
    plt.axhline(y=1, color='#4daf4a', linestyle='--')

    plt.legend()
    plt.xlabel('Metabolites')
    plt.xticks([])
    plt.ylabel('Correlation with verified integration')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)


if __name__ == '__main__':
    """plot_MetCohort_cor(159, 162, 'C5H13NO', 'C[13]C4H13NO')
    # plot_xcms_cor(329, 346, 'C5H13NO', 'C[13]C4H13NO')

    plot_MetCohort_cor(194, 199, 'C5H11NO2', 'C[13]C4H11NO2')
    # plot_xcms_cor(440, 443, 'C5H11NO2', 'C[13]C4H11NO2')

    plot_MetCohort_cor(263, 270, 'C9H7N', 'C[13]C8H7N')
    # plot_xcms_cor(539, 558, 'C9H7N', 'C[13]C8H7N')

    plot_MetCohort_cor(267, 274, 'C8H19N', 'C[13]C7H19N')
    # plot_xcms_cor(545, 562, 'C8H19N', 'C[13]C7H19N')

    plot_MetCohort_cor(511, 521, 'C10H9NO2', 'C[13]C9H9NO2')
    # plot_xcms_cor(972, 987, 'C10H9NO2', 'C[13]C9H9NO2')

    plot_MetCohort_cor(783, 791, 'C11H11NO3', 'C[13]C10H11NO3')
    # plot_xcms_cor(1329, 1337, 'C11H11NO3', 'C[13]C10H11NO3')

    plot_MetCohort_cor(842, 853, 'C15H17N', 'C[13]C14H17N')
    # plot_xcms_cor(1393, 1402, 'C15H17N', 'C[13]C14H17N')

    plot_MetCohort_cor(1033, 1041, 'C2H2KNaO5S2', 'C[13]CH2KNaO5S2')
    # plot_xcms_cor(1599, 1609, 'C2H2KNaO5S2', 'C[13]CH2KNaO5S2')

    plot_MetCohort_cor(1096, 1100, 'C10H22O6', '[13]CC9H22O6')
    # plot_xcms_cor(1657, 1660, 'C10H22O6', '[13]CC9H22O6')

    plot_MetCohort_cor(1192, 1203, 'C14H18N2O2', '[13]CC13H18N2O2')
    # plot_xcms_cor(1752, 1764, 'C14H18N2O2', '[13]CC13H18N2O2')

    plot_MetCohort_cor(1354, 1366, 'C13H25NO4', '[13]CC12H25NO4')
    # plot_xcms_cor(1895, 1912, 'C13H25NO4', '[13]CC12H25NO4')

    plot_MetCohort_cor(183, 300, 'C4H7N3O', 'C4H6N3O+Na')
    # plot_xcms_cor(417, 602, 'C4H7N3O', 'C4H6N3O+Na')

    plot_MetCohort_cor(388, 533, 'C6H9N3O2', 'C6H8N3O2+Na')
    # plot_xcms_cor(767, 996, 'C6H9N3O2', 'C6H8N3O2+Na')

    plot_MetCohort_cor(773, 979, 'C11H12N2O2', 'C11H11N2O2+Na')
    # plot_xcms_cor(1323, 1540, 'C11H12N2O2', 'C11H11N2O2+Na')

    plot_MetCohort_cor(1063, 1320, 'C13H21N3O', 'C13H20N3O+Na')
    # plot_xcms_cor(1635, 1859, 'C13H21N3O', 'C13H20N3O+Na')

    plot_MetCohort_cor(1241, 1512, 'C17H18N2', 'C17H17N2+Na')
    # plot_xcms_cor(1790, 2042, 'C17H18N2', 'C17H17N2+Na')

    plot_MetCohort_cor(1252, 1520, 'C17H17NO', 'C17H16NO+Na')
    # plot_xcms_cor(1801, 2052, 'C17H17NO', 'C17H16NO+Na')

    plot_MetCohort_cor(1413, 1694, 'C13H16N2O4', 'C13H15N2O4+Na')
    # plot_xcms_cor(1950, 2209, 'C13H16N2O4', 'C13H15N2O4+Na')

    plot_MetCohort_cor(2157, 2470, 'C16H18N2O5', 'C16H18N2O5+Na')
    # plot_xcms_cor(2613, 2882, 'C16H18N2O5', 'C16H18N2O5+Na')

    plot_MetCohort_cor(265, 341, 'C6H11NO2', 'C6H11NO2+NH3')
    # plot_xcms_cor(540, 689, 'C6H11NO2', 'C6H11NO2+NH3')

    plot_MetCohort_cor(611, 774, 'C11H9NO2', 'C11H9NO2+NH3')
    # plot_xcms_cor(1082, 1323, 'C11H9NO2', 'C11H9NO2+NH3')

    plot_MetCohort_cor(838, 993, 'C13H22O2', 'C13H22O2+NH3')
    # plot_xcms_cor(1389, 1563, 'C13H22O2', 'C13H22O2+NH3')

    plot_MetCohort_cor(1308, 1524, 'C14H24O2S', 'C14H24O2S+NH3')
    # plot_xcms_cor(1849, 2055, 'C14H24O2S', 'C14H24O2S+NH3')

    plot_MetCohort_cor(1649, 1873, 'C12H26O7', 'C12H26O7+NH3')
    # plot_xcms_cor(2167, 2362, 'C12H26O7', 'C12H26O7+NH3')

    plot_MetCohort_cor(1904, 2157, 'C16H15NO5', 'C16H15NO5+NH3')
    # plot_xcms_cor(2385, 2613, 'C16H15NO5', 'C16H15NO5+NH3')

    plot_MetCohort_cor(2471, 2706, 'C14H28O9', 'C14H28O9+NH3')
    # plot_xcms_cor(2884, 3086, 'C14H28O9', 'C14H28O9+NH3')

    plot_MetCohort_cor(2472, 2708, 'C22H28O3', 'C22H28O3+NH3')
    # plot_xcms_cor(2887, 3087, 'C22H28O3', 'C22H28O3+NH3')"""

    # compare_integration("D:\\Experiments\\peakmat_evaluation\\XG\\xgym_pos_manual.xlsx")

    plot_correlation(os.path.join(picture_dir, 'quantification_correlation.xlsx'),
                     os.path.join(picture_dir, 'integration_correlation_renamed.tiff'))

