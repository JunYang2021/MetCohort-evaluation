# Figure 4d, Figure 4g, Figure s9,
# Figure s10, Figure s11b, Figure s12,
# Figure s13, Figure s14, Figure s15

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from pyopenms import *
import math
from scipy.stats import gaussian_kde
from BTrees.OOBTree import OOBTree

font = {'family': 'Arial', 'size': 7}
plt.rc('font', **font)


def compare_features(file_path, delta_mz=0.02, delta_rt=5):
    true_features = pd.read_csv('mix39 true features.csv')
    comparison_file = pd.read_csv(file_path)
    comparison_file_name = os.path.basename(file_path)

    true_features[comparison_file_name] = 0

    for index, row in true_features.iterrows():
        mz_value = row['mz']
        rt_value = row['rt']

        # Check if there's a match within the delta values in comparison_file
        condition = (
                (comparison_file['mz'] >= mz_value - delta_mz) &
                (comparison_file['mz'] <= mz_value + delta_mz) &
                (comparison_file['rt'] >= rt_value - delta_rt) &
                (comparison_file['rt'] <= rt_value + delta_rt)
        )
        if any(condition):
            true_features.at[index, comparison_file_name] = 1

    # Save the modified 'true features.csv'
    true_features.to_csv('mix39 true features.csv', index=False)


def compare_feature_excel(file_path, delta_mz=0.02, delta_rt=5):
    excel_data = pd.read_excel('feature detection evaluation MIX39.xlsx',
                               sheet_name=0)
    csv_data = pd.read_csv(file_path)

    csv_data['name'] = ''

    for index, row in excel_data.iterrows():
        mz_value = row['mz']
        rt_value = row['rt']

        condition = (
                (csv_data['mz'] >= mz_value - delta_mz) &
                (csv_data['mz'] <= mz_value + delta_mz) &
                (csv_data['rt'] >= rt_value - delta_rt) &
                (csv_data['rt'] <= rt_value + delta_rt)
        )

        matching_rows = csv_data[condition]

        if not matching_rows.empty:
            matching_rows['distance'] = (matching_rows['rt'] - rt_value).abs()

            matching_rows = matching_rows.sort_values(by='distance')

            # iterate over the sorted rows and find the first row without a name
            for _, match_row in matching_rows.iterrows():
                if match_row['name']:
                    print(row)
                if not match_row['name']:
                    csv_data.at[match_row.name, 'name'] = row['name']
                    break

    with pd.ExcelWriter('feature detection evaluation MIX39.xlsx',
                        engine='openpyxl',
                        mode='a') as writer:
        # only write mz, rt, and name columns to the Excel file
        csv_data[['mz', 'rt', 'name']].to_excel(writer,
                                                sheet_name=os.path.basename(file_path).split('.')[0],
                                                index=False)


def plot_xic(file_list, color_list, mz_range, time_range, picture_name=None):
    # fig, ax = plt.subplots(figsize=(5.8, 2.3))
    fig, ax = plt.subplots(figsize=(4, 2))

    plt.rcParams['font.sans-serif'] = "Arial"
    # plt.rcParams['font.size'] = 12
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = "sans-serif"

    for file, color in zip(file_list, color_list):
        exp = MSExperiment()
        MzMLFile().load(file, exp)
        time = []
        intensity = []
        for scan in exp:
            if scan.getMSLevel() == 1:
                if time_range[0] < scan.getRT() < time_range[1]:
                    mz, i = scan.get_peaks()
                    mask = (mz > mz_range[0]) & (mz < mz_range[1])
                    time.append(scan.getRT())
                    intensity.append(i[mask].sum())
                elif scan.getRT() > time_range[1]:
                    break
        plt.plot(time, intensity, color=color, linewidth=0.7)
        # plt.plot(time, intensity, color=color, linewidth=1.7)

    # plt.ticklabel_format(axis='both', style='sci', scilimits=(4, 4))
    plt.xlabel('Retention Time (s)')
    plt.ylabel('Intensity')
    plt.title('m/z: {0}-{1}'.format(mz_range[0], mz_range[1]))
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    plt.locator_params(axis='x', nbins=4)
    plt.locator_params(axis='y', nbins=4)

    ax.set_xlim(time_range)
    # plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tight_layout()

    if picture_name is not None:
        plt.savefig(picture_name, dpi=300)
    # plt.show()


def compute_ei(i_array):
    max_intensity = max(i_array)
    max_index = i_array.index(max_intensity)

    variant_intensities = []
    # Left side of the peak
    for i in range(1, max_index):
        if i_array[i] < i_array[i - 1]:
            variant_intensities.append(i_array[i - 1] - i_array[i])
    # Right side of the peak
    for i in range(max_index + 1, len(i_array)):
        if i_array[i] > i_array[i - 1]:
            variant_intensities.append(i_array[i] - i_array[i - 1])

    if not variant_intensities:
        entropy_index = 0
    else:
        total = max_intensity + sum(variant_intensities)
        p_max = max_intensity / total
        p_variants = [vi / total for vi in variant_intensities]
        entropy_index = -p_max * math.log(p_max) - sum([pv * math.log(pv) for pv in p_variants])

    return entropy_index


def compute_ei_height(msexperiment, feature_mz, feature_rt, delta_mz=0.01, delta_rt=5):
    intensity = []
    for scan in msexperiment:
        if scan.getMSLevel() == 1:
            if feature_rt - delta_rt < scan.getRT() < feature_rt + delta_rt:
                mz, i = scan.get_peaks()
                mask = (mz > feature_mz - delta_mz) & (mz < feature_mz + delta_mz)
                intensity.append(i[mask].sum())
            elif scan.getRT() > feature_rt + delta_rt:
                break

    max_intensity = max(intensity)
    entropy_index = compute_ei(intensity)

    return entropy_index, max_intensity


def plot_xic_ei(file, mz_range, time_range, picture_name):
    fig, ax = plt.subplots(figsize=(3, 2))
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    time = []
    intensity = []
    for scan in exp:
        if scan.getMSLevel() == 1:
            if time_range[0] < scan.getRT() < time_range[1]:
                mz, i = scan.get_peaks()
                mask = (mz > mz_range[0]) & (mz < mz_range[1])
                time.append(scan.getRT())
                intensity.append(i[mask].sum())
            elif scan.getRT() > time_range[1]:
                break

    ei = compute_ei(intensity)
    plt.plot(time, intensity, linewidth=0.7)
    plt.xlabel('Retention Time (s)')
    plt.ylabel('Intensity')
    plt.title('Entropy index: {:.2f}'.format(ei))
    plt.locator_params(axis='x', nbins=4)
    plt.locator_params(axis='y', nbins=4)

    ax.set_xlim(time_range)
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tight_layout()
    plt.savefig(picture_name, dpi=300)


class UniqueFeature:
    def __init__(self, mz, rt, file, id):
        self.mz = [mz]
        self.rt = [rt]
        self.file = [file]
        self.meanmz = mz
        self.meanrt = rt
        self.unique_file = set(self.file)
        self.id = [id]


def read_excel_features(file):
    try:
        df = pd.read_excel(file, engine="openpyxl")
        data = [(row["m/z"], row["rt"], int(row['id'])) for index, row in df.iterrows()]
        print("Reading {0} features from {1}".format(len(data), file))
        return data
    except FileNotFoundError:
        print(f"The file {file} does not exist. Please check the filename and try again.")
    except pd.errors.ParserError:
        print(f"There was an error parsing the file {file}. Please check the file format and try again.")
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")


def merge_features(files_list, delta_mz=0.01, delta_rt=3):
    """

    :param files_list:
    :param delta_mz: m/z deviation allowed
    :param delta_rt: rt deviation allowed(s)
    :return:
    """
    file_names = [os.path.splitext(os.path.basename(i))[0] for i in files_list]

    unique_feature_tree = OOBTree()
    for i, file in enumerate(files_list):
        feature_set = read_excel_features(file)
        if i == 0:
            unique_feature_tree[feature_set[0][0]] = UniqueFeature(feature_set[0][0],
                                                                   feature_set[0][1],
                                                                   file_names[i],
                                                                   feature_set[0][2])
        for j, feature in enumerate(feature_set):
            if i == 0 and j == 0:
                continue
            append_or_not = False
            for inrange_mz, inrange_item in unique_feature_tree.items(min=feature[0] - delta_mz,
                                                                      max=feature[0] + delta_mz):
                if abs(inrange_item.meanrt - feature[1]) < delta_rt:
                    inrange_item.mz.append(feature[0])
                    inrange_item.rt.append(feature[1])
                    inrange_item.id.append(feature[2])
                    inrange_item.file.append(file_names[i])
                    inrange_item.meanmz = sum(inrange_item.mz) / len(inrange_item.mz)
                    inrange_item.meanrt = sum(inrange_item.rt) / len(inrange_item.rt)
                    inrange_item.unique_file = set(inrange_item.file)
                    append_or_not = True
                    break
            if not append_or_not:
                if feature[0] in unique_feature_tree:
                    unique_feature_tree[feature[0] + 0.0000001] = UniqueFeature(feature[0], feature[1], file_names[i],
                                                                                feature[2])
                else:
                    unique_feature_tree[feature[0]] = UniqueFeature(feature[0], feature[1], file_names[i], feature[2])
    return unique_feature_tree


def print_feature_tree(feature_tree, output_file):
    data = []

    for mz, unique_feature in feature_tree.items():
        data.append({
            "mz": unique_feature.meanmz,
            "rt": unique_feature.meanrt,
            "mz_list": str(unique_feature.mz),
            "rt_list": str(unique_feature.rt),
            "file_list": str(unique_feature.file),
            "unique_file": str(unique_feature.unique_file),
            "id": str(unique_feature.id)
        })

    df = pd.DataFrame(data)
    df.to_excel(output_file, index=False, engine='openpyxl')


# From asari: mSelectivity (https://github.com/shuzhao-li-lab/asari/blob/main/asari/mass_functions.py)
def calculate_selectivity(sorted_mz_list, std_ppm=5):
    '''
    To calculate m or d-selectivity for a list of m/z values,
    which can be all features in an experiment or a database.
    The mass selectivity between two m/z values is defined as:
    (1 - Probability(confusing two peaks)), further formalized as an exponential model:
    P = exp( -x/std_ppm ),
    whereas x is ppm distance between two peaks,
    std_ppm standard deviation of ppm between true peaks and theoretical values, default at 5 pmm.

    The selectivity value is between (0, 1), close to 1 meaning high selectivity.
    If multiple adjacent peaks are present, we multiply the selectivity scores.
    It is good approximation by considering 2 lower and 2 higher neighbors here.

    Parameters
    ----------
    sorted_mz_list: list
        a list of m/z values, sorted from low to high, length > 3.
    std_ppm: float, optional, default: 5
        mass resolution in ppm (part per million).

    Returns
    -------
    A list of selectivity values, in matched order as the input m/z list.


    Note
    ----
    ppm is actually dependent on m/z, not an ideal method.
    But it's in common practice and good enough approximation.
    '''

    def __sel__(x, std_ppm=std_ppm):
        if x > 100:  # too high, not bother
            return 1
        elif x < 0.1:
            return 0
        else:
            return 1 - np.exp(-x / std_ppm)

    mz_list = np.array(sorted_mz_list)
    ppm_distances = 1000000 * (mz_list[1:] - mz_list[:-1]) / mz_list[:-1]
    # first two MassTraces
    selectivities = [
        __sel__(ppm_distances[0]) * __sel__(ppm_distances[0] + ppm_distances[1]),
        __sel__(ppm_distances[0]) * __sel__(ppm_distances[1]
                                            ) * __sel__(ppm_distances[1] + ppm_distances[2]),
    ]
    for ii in range(2, mz_list.size - 2):
        selectivities.append(
            __sel__(ppm_distances[ii - 2] + ppm_distances[ii - 1]) * __sel__(
                ppm_distances[ii - 1]) * __sel__(ppm_distances[ii]) * __sel__(
                ppm_distances[ii] + ppm_distances[ii + 1])
        )
    # last two MassTraces
    selectivities += [
        __sel__(ppm_distances[-3] + ppm_distances[-2]) * __sel__(
            ppm_distances[-2]) * __sel__(ppm_distances[-1]),
        __sel__(ppm_distances[-2] + ppm_distances[-1]) * __sel__(ppm_distances[-1]),
    ]

    return selectivities


def compute_ei_height_excel(excel_path, qc_paths, delta_mz=0.01, delta_rt=5):
    df = pd.read_excel(excel_path, engine='openpyxl')

    all_entropy_indices = []
    all_heights = []

    for qc_path in qc_paths:
        exp = MSExperiment()
        MzMLFile().load(qc_path, exp)

        entropy_indices = []
        heights = []

        for _, row in df.iterrows():
            feature_mz = row['mz']
            feature_rt = row['rt']
            entropy_index, height = compute_ei_height(exp, feature_mz, feature_rt, delta_mz, delta_rt)
            entropy_indices.append(entropy_index)
            heights.append(height)

        all_entropy_indices.append(entropy_indices)
        all_heights.append(heights)

    # Compute the median of entropy index and height across all files
    median_entropy_indices = [np.median([ei[i] for ei in all_entropy_indices]) for i in range(len(df))]
    median_heights = [np.median([h[i] for h in all_heights]) for i in range(len(df))]

    df['entropy index'] = median_entropy_indices
    df['height'] = median_heights

    df.to_excel(excel_path, index=False, engine='openpyxl')


def plot_ei_height(excel_path, labels, out_path=None):
    plt.figure(figsize=(5, 4))

    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.size'] = 12
    plt.rcParams['font.family'] = "sans-serif"

    df = pd.read_excel(excel_path, engine='openpyxl')
    df['height_adjusted'] = df['height'].apply(lambda x: 0.1 if x <= 0 else x)

    unique_classes = df['unique_file'].unique()

    unique_classes[2], unique_classes[1] = unique_classes[1], unique_classes[2]
    # unique_classes[0], unique_classes[1], unique_classes[2] = unique_classes[1], unique_classes[0], unique_classes[2]
    print(unique_classes)
    colors = ["#e41a1c", "#377eb8", "#4daf4a"]
    color_map = {class_val: colors[i] for i, class_val in enumerate(unique_classes)}
    legend_map = {unique_classes[i]: label for i, label in enumerate(labels)}

    for class_val, color in color_map.items():
        subset = df[df['unique_file'] == class_val]
        plt.scatter(subset['entropy index'], subset['height_adjusted'], c=color, label=legend_map[class_val], s=5,
                    alpha=0.6)

    plt.xlabel('Entropy index')
    plt.ylabel('Height')
    plt.yscale('log')
    plt.ylim(30, 10 ** 7.7)
    plt.yticks([10 ** i for i in range(2, 8)])
    plt.gca().tick_params(axis='y', which='minor', size=0)

    plt.legend()
    plt.tight_layout()
    if out_path:
        plt.savefig(out_path, dpi=300)


def plot_xic_ei_example(file, _mz, _rt, picture_name):
    fig, ax = plt.subplots(figsize=(3, 2))
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    plot_time = []      # span 20 s
    plot_intensity = []
    cal_time = []       # span 10 s
    cal_intensity = []
    for scan in exp:
        if scan.getMSLevel() == 1:
            if _rt - 10 < scan.getRT() < _rt + 10:
                mz, i = scan.get_peaks()
                mask = (mz > _mz - 0.01) & (mz < _mz + 0.01)
                plot_time.append(scan.getRT())
                plot_intensity.append(i[mask].sum())
                if _rt - 5 < scan.getRT() < _rt + 5:
                    cal_time.append(scan.getRT())
                    cal_intensity.append(i[mask].sum())
            elif scan.getRT() > _rt + 10:
                break

    ei = compute_ei(cal_intensity)
    plt.plot(plot_time, plot_intensity, linewidth=1.1)
    plt.xlabel('Retention Time (s)')
    plt.ylabel('Intensity')
    plt.title('Entropy index: {:.2f}'.format(ei))
    plt.locator_params(axis='x', nbins=4)
    plt.locator_params(axis='y', nbins=4)

    ax.set_xlim((plot_time[0], plot_time[-1]))
    plt.tick_params(axis='both', which='major', labelsize=6)
    plt.tight_layout()
    plt.savefig(picture_name, dpi=300)


def plot_mselectivity(peak_list_file, msrange=None, out_putfile=None):
    df = pd.read_excel(peak_list_file, engine="openpyxl")
    mz_list = df['m/z'].sort_values().tolist()
    ms = calculate_selectivity(mz_list)
    density = gaussian_kde(ms)
    if not msrange:
        x_range = np.linspace(min(ms), max(ms), 500)
    else:
        x_range = np.linspace(msrange[0], msrange[1], 500)
    plt.figure(figsize=(4, 3))
    plt.plot(x_range, density(x_range))
    plt.xlabel('mSelectivity')
    plt.ylabel('Kernel density')
    plt.savefig(out_putfile, dpi=300)


def plot_xic_exp(exp_list, color_list, feature_mz, feature_rt, picture_name=None):
    mz_range = (feature_mz - 0.01, feature_mz + 0.01)
    time_range = (feature_rt - 5, feature_rt + 5)

    # fig, ax = plt.subplots(figsize=(5.8, 2.3))
    fig, ax = plt.subplots(figsize=(1.5, 1))

    plt.rcParams['font.sans-serif'] = "Arial"
    # plt.rcParams['font.size'] = 12
    plt.rcParams['font.size'] = 4
    plt.rcParams['font.family'] = "sans-serif"

    ei, height = [], []
    for exp in exp_list:
        _ei, _height = compute_ei_height(exp, feature_mz, feature_rt, delta_rt=10)
        ei.append(_ei)
        height.append(_height)

    for exp, color in zip(exp_list, color_list):
        time = []
        intensity = []
        for scan in exp:
            if scan.getMSLevel() == 1:
                if time_range[0] < scan.getRT() < time_range[1]:
                    mz, i = scan.get_peaks()
                    mask = (mz > mz_range[0]) & (mz < mz_range[1])
                    time.append(scan.getRT())
                    intensity.append(i[mask].sum())
                elif scan.getRT() > time_range[1]:
                    break
        plt.plot(time, intensity, color=color, linewidth=0.7)
        # plt.plot(time, intensity, color=color, linewidth=1.7)

    # plt.ticklabel_format(axis='both', style='sci', scilimits=(4, 4))
    plt.xlabel('Retention Time (s)')
    plt.ylabel('Intensity')
    plt.title('m/z: {:.4f}-{:.4f}'.format(mz_range[0], mz_range[1]))
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    plt.locator_params(axis='x', nbins=3)
    plt.locator_params(axis='y', nbins=3)

    ax.set_xlim(time_range)
    # plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='major', labelsize=4)

    x_pos = ax.get_xlim()[0] + 0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0])  # 5% from the left boundary
    y_pos_max = ax.get_ylim()[1]
    y_pos_step = 0.15 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.text(x_pos, y_pos_max - y_pos_step, f"Entropy index = {np.median(ei):.2f}", fontsize=3, color='black')
    ax.text(x_pos, y_pos_max - 2 * y_pos_step, f"Feature height = {np.median(height):.2f}", fontsize=3, color='black')
    plt.tight_layout()

    if picture_name is not None:
        plt.savefig(picture_name, dpi=300)
    plt.close()


if __name__ == '__main__':
    """sort_mz = np.array([100, 101, 105, 108.5, 108.6, 120, 120.01])
    mselectivities = calculate_selectivity(sort_mz)
    print(mselectivities)"""

    """pd.set_option('mode.chained_assignment', None)
    compare_feature_excel(file_path='../gpw_mix39/xcms-16253.csv', delta_mz=0.02, delta_rt=10)
    print('*' * 50)
    compare_feature_excel(file_path='../gpw_mix39/mzmine-33721.csv', delta_mz=0.02, delta_rt=10)
    print('*' * 50)
    compare_feature_excel(file_path='../gpw_mix39/msdial-20117.csv', delta_mz=0.02, delta_rt=10)
    print('*' * 50)
    compare_feature_excel(file_path='../gpw_mix39/peakmat_v3_15910.csv', delta_mz=0.02, delta_rt=10)  # Untargeted

    compare_features(file_path='../gpw_mix39/xcms-16253.csv', delta_mz=0.02, delta_rt=10)
    compare_features(file_path='../gpw_mix39/mzmine-33721.csv', delta_mz=0.02, delta_rt=10)
    compare_features(file_path='../gpw_mix39/msdial-20117.csv', delta_mz=0.02, delta_rt=10)
    compare_features(file_path='../gpw_mix39/peakmat_v3_15910.csv', delta_mz=0.02, delta_rt=10)"""

    # Plot of MIX39
    """mix39 = 'F:\\Test_Peakmat\\gpw_mix39\\mzml_data'
    mix39_names = [os.path.join(mix39, f) for f in os.listdir(mix39) if f.endswith('.mzML')]
    mix39_names = [os.path.normpath(i) for i in mix39_names]
    color_list = plt.cm.turbo(np.linspace(0, 1, 20))
    cpd = [dict(name="4,4'-DiaminodiphenylEther", mzrange=(129.0022, 129.0222), timerange=(36.8, 56.8)),
           dict(name='Secbumeton', mzrange=(226.1568, 226.1768), timerange=(162, 182)),
           dict(name='Clenbuterol', mzrange=(277.0774, 277.0974), timerange=(152, 172)),
           dict(name='Difenoxuron', mzrange=(287.1296, 287.1496), timerange=(252, 272))]

    for cd in cpd:
        plot_xic(mix39_names,
                 color_list,
                 cd['mzrange'],
                 cd['timerange'],
                 r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\{0}.tiff'.format(
                     cd['name']))"""

    # Example of entropy index
    ep_file = r'G:\Test_Peakmat\MSV000084790\mzml\20180126_7_3-1A_1_neg_100.mzML'
    plot_xic_ei_example(ep_file, 85.0516, 356,
                        r'ei_low_example.tiff')
    plot_xic_ei_example(ep_file, 250.804, 48,
                        r'ei_middle_example.tiff')
    plot_xic_ei_example(ep_file, 269.2324, 306,
                        r'ei_middle_2_example.tiff')
    plot_xic_ei_example(ep_file, 183.5352, 351,
                        r'ei_high_example.tiff')

    # MSV000084790 data false positive evaluation
    """# peakmat vs xcms
    file1 = 'F:\\Test_Peakmat\\MSV000084790\\PEAKMAT-RES\\peakmat-v3-7066.xlsx'
    file2 = 'F:\\Test_Peakmat\\MSV000084790\\XCMS-RES\\xcms-4095.xlsx'

    outfile = 'F:\\Test_Peakmat\\MSV000084790\\peakmat_vs_xcms_20240523.xlsx'

    feature_tree = merge_features([file1, file2], delta_mz=0.01, delta_rt=5)
    area1, area2, n12 = [0 for i in range(3)]
    for _, feature in feature_tree.items():
        if 'peakmat-v3-7066' in feature.unique_file:
            area1 += 1
        if 'xcms-4095' in feature.unique_file:
            area2 += 1
        if {'peakmat-v3-7066', 'xcms-4095'}.issubset(feature.unique_file):
            n12 += 1
    print(area1, area2, n12)
    print_feature_tree(feature_tree, outfile)

    # peakmat vs mzmine
    file1 = 'F:\\Test_Peakmat\\MSV000084790\\PEAKMAT-RES\\peakmat-v3-7066.xlsx'
    file2 = 'F:\\Test_Peakmat\\MSV000084790\\MZMINE-RES\\mzmine-7699.xlsx'

    outfile = 'F:\\Test_Peakmat\\MSV000084790\\peakmat_vs_mzmine_20240523.xlsx'

    feature_tree = merge_features([file1, file2], delta_mz=0.01, delta_rt=5)
    area1, area2, n12 = [0 for i in range(3)]
    for _, feature in feature_tree.items():
        if 'peakmat-v3-7066' in feature.unique_file:
            area1 += 1
        if 'mzmine-7699' in feature.unique_file:
            area2 += 1
        if {'peakmat-v3-7066', 'mzmine-7699'}.issubset(feature.unique_file):
            n12 += 1
    print(area1, area2, n12)
    print_feature_tree(feature_tree, outfile)

    # peakmat vs asari
    file1 = 'F:\\Test_Peakmat\\MSV000084790\\PEAKMAT-RES\\peakmat-v3-7066.xlsx'
    file2 = 'F:\\Test_Peakmat\\MSV000084790\\asari\\output_asari_project_523153855\\asari_msv84790-1070.xlsx'

    outfile = 'F:\\Test_Peakmat\\MSV000084790\\peakmat_vs_asari_20240523.xlsx'

    feature_tree = merge_features([file1, file2], delta_mz=0.01, delta_rt=5)
    area1, area2, n12 = [0 for i in range(3)]
    for _, feature in feature_tree.items():
        if 'peakmat-v3-7066' in feature.unique_file:
            area1 += 1
        if 'asari_msv84790-1070' in feature.unique_file:
            area2 += 1
        if {'peakmat-v3-7066', 'asari_msv84790-1070'}.issubset(feature.unique_file):
            n12 += 1
    print(area1, area2, n12)
    print_feature_tree(feature_tree, outfile)"""

    """excel_path1 = 'F:\\Test_Peakmat\\MSV000084790\\peakmat_vs_mzmine_20240523.xlsx'
    excel_path2 = 'F:\\Test_Peakmat\\MSV000084790\\peakmat_vs_xcms_20240523.xlsx'
    excel_path3 = 'F:\\Test_Peakmat\\MSV000084790\\peakmat_vs_asari_20240523.xlsx'
    qc_path = ['F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_QC-EI_neg_112.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_14_5-1_2_neg_065.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_117_87-3A_4_neg_228.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_121_91-3A_1_neg_085.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_160_119-1B_1_neg_132.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_199_141-3A_1_neg_104.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_201_141-3A_3_neg_098.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_QC-EI_neg_068.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_QC-EI_neg_189.mzML',
               'F:\\Test_Peakmat\\MSV000084790\\mzml\\20180126_QC-EI_neg_245.mzML']
    compute_ei_height_excel(excel_path1, qc_path, delta_rt=5)
    compute_ei_height_excel(excel_path2, qc_path, delta_rt=5)
    compute_ei_height_excel(excel_path3, qc_path, delta_rt=5)"""

    """excel_path1 = 'D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_mzmine_20240523.xlsx'
    excel_path2 = 'D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_xcms_20240523.xlsx'
    excel_path3 = 'D:\\Experiments\\peakmat_evaluation\\msv84790\\peakmat_vs_asari_20240523.xlsx'
    labels1 = ['MetCohort & MZmine 3', 'Only MZmine 3', 'Only MetCohort']
    labels2 = ['MetCohort & XCMS', 'Only XCMS', 'Only MetCohort']
    labels3 = ['MetCohort & asari', 'Only asari', 'Only MetCohort']
    # plot_ei_height(excel_path1, labels1, out_path=r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\MSV84790_metcohort_mzmine.tiff')
    plot_ei_height(excel_path2, labels2, out_path=r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\MSV84790_metcohort_xcms.tiff')
    # plot_ei_height(excel_path3, labels3,
    #                out_path=r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\MSV84790_metcohort_asari.tiff')"""

    """peakmat_list = r'D:/Experiments/peakmat_evaluation/msv84790/peakmat-v3-7066.xlsx'
    xcms_list = r'D:/Experiments/peakmat_evaluation/msv84790/xcms-4095.xlsx'
    mzmine_list = r'D:/Experiments/peakmat_evaluation/msv84790/mzmine-7699.xlsx'
    asari_list = r'D:/Experiments/peakmat_evaluation/msv84790/asari_msv84790-1070.xlsx'
    plot_mselectivity(peakmat_list, out_putfile=r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\MSV84790_peakmat_ms.tiff')
    plot_mselectivity(xcms_list, out_putfile=r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\MSV84790_xcms_ms.tiff')
    plot_mselectivity(mzmine_list,
                      out_putfile=r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\MSV84790_mzmine_ms.tiff')
    plot_mselectivity(asari_list,
                      out_putfile=r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\MSV84790_asari_ms.tiff')"""

    ## plot of MSV000084790  features detected by different tools
    """msv90 = 'F:\\Test_Peakmat\\MSV000084790\\mzml'
    msv90_names = [os.path.join(msv90, f) for f in os.listdir(msv90) if f.endswith('.mzML')][-20:]
    msv90_names = [os.path.normpath(i) for i in msv90_names]
    exps = []
    for file in msv90_names:
        exp = MSExperiment()
        MzMLFile().load(file, exp)
        exps.append(exp)

    color_list = plt.cm.turbo(np.linspace(0, 1, 20))
    cpd = [dict(name='overlap_p_x_1', featuremz=223.096699921324, featurert=549.26),
           dict(name="overlap_p_x_2", featuremz=341.225514651028, featurert=434.24),
           dict(name='overlap_p_x_3', featuremz=424.77926343749, featurert=369.67),
           dict(name='overlap_p_x_4', featuremz=560.842920646819, featurert=391.70),
           dict(name='x_only_p_1', featuremz=104.941741943359, featurert=459.36),
           dict(name='x_only_p_2', featuremz=133.085723876953, featurert=356.68),
           dict(name='x_only_p_3', featuremz=342.1787415, featurert=432.86),
           dict(name='x_only_p_4', featuremz=459.011657714844, featurert=325.01),
           dict(name='x_only_x_1', featuremz=87.0443502254364, featurert=484.29),
           dict(name='x_only_x_2', featuremz=338.7035924, featurert=392.48),
           dict(name='x_only_x_3', featuremz=443.3348500002, featurert=761.29),
           dict(name='x_only_x_4', featuremz=649.451841174876, featurert=721.63),
           dict(name='overlap_p_m_1', featuremz=87.2759017393359, featurert=365.66),
           dict(name="overlap_p_m_2", featuremz=224.076032193848, featurert=99.85),
           dict(name='overlap_p_m_3', featuremz=361.198214903125, featurert=661.23),
           dict(name='overlap_p_m_4', featuremz=485.359711317578, featurert=268.95),
           dict(name='m_only_p_1', featuremz=91.0584106445312, featurert=46.93),
           dict(name='m_only_p_2', featuremz=294.241851806641, featurert=487.64),
           dict(name='m_only_p_3', featuremz=402.731323242187, featurert=393.91),
           dict(name='m_only_p_4', featuremz=466.241943359375, featurert=396.64),
           dict(name='m_only_m_1', featuremz=90.05532375, featurert=104.07),
           dict(name='m_only_m_2', featuremz=217.9599252, featurert=131.21),
           dict(name='m_only_m_3', featuremz=372.0641774, featurert=293.70),
           dict(name='m_only_m_4', featuremz=512.2861684, featurert=331.13),
           dict(name='overlap_p_a_1', featuremz=114.819008569336, featurert=324.953),
           dict(name="overlap_p_a_2", featuremz=254.824259375, featurert=49.45),
           dict(name='overlap_p_a_3', featuremz=373.714544628906, featurert=361.13),
           dict(name='overlap_p_a_4', featuremz=497.235660522461, featurert=324.953),
           dict(name='a_only_p_1', featuremz=86.0561828613281, featurert=292.75),
           dict(name='a_only_p_2', featuremz=224.864654541016, featurert=46.25),
           dict(name='a_only_p_3', featuremz=370.122039794922, featurert=275.20),
           dict(name='a_only_p_4', featuremz=498.313385009766, featurert=447.90),
           dict(name='a_only_a_1', featuremz=111.0204, featurert=45.49),
           dict(name='a_only_a_2', featuremz=341.2128, featurert=340.97),
           dict(name='a_only_a_3', featuremz=382.9832, featurert=326.64),
           dict(name='a_only_a_4', featuremz=564.3591, featurert=325.76)
           ]

    for cd in cpd:
        plot_xic_exp(exps,
                     color_list,
                     cd['featuremz'],
                     cd['featurert'],
                     r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure4\msv790\{0}.tiff'.format(cd['name']))"""

    """# For SP61 data evaluation
    # peakmat vs xcms
    file1 = 'D:\\Experiments\\peakmat_evaluation\\SP61\\peakmat_v3_5358.xlsx'
    file2 = 'D:\\Experiments\\peakmat_evaluation\\SP61\\xcms-6807.xlsx'

    outfile = 'D:\\Experiments\\peakmat_evaluation\\SP61\\peakmat_vs_xcms_20240525.xlsx'

    feature_tree = merge_features([file1, file2], delta_mz=0.01, delta_rt=5)
    area1, area2, n12 = [0 for i in range(3)]
    for _, feature in feature_tree.items():
        if 'peakmat_v3_5358' in feature.unique_file:
            area1 += 1
        if 'xcms-6807' in feature.unique_file:
            area2 += 1
        if {'peakmat_v3_5358', 'xcms-6807'}.issubset(feature.unique_file):
            n12 += 1
    print(area1, area2, n12)        # 5172 6807 3584
    print_feature_tree(feature_tree, outfile)"""