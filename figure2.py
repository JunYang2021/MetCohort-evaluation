from src.file_processing import MsFile, ms1scan
from src.dtw import dtw, loess_dtw
import os
import matplotlib.pyplot as plt
import numpy as np
from src.get_chromatograms import roa_construction, get_eic_array
from src.local_match import local_match
from scipy import signal
from scipy.stats import t

font = {'family': 'Arial', 'size': 7}
plt.rc('font', **font)


# perform dtw on a given file list (generate delta rt plot and aligned files)
def alignment_only_dtw(input_files, new_dir, figure_file):
    msfiles = []
    for file in input_files:
        msfiles.append(MsFile(file))

    id_ref = 0
    ref_file = input_files[id_ref]
    fig, ax = plt.subplots(figsize=(3, 2.5))
    colors = plt.cm.rainbow(np.linspace(0, 1, len(input_files)))
    for _i, file in enumerate(input_files):
        if file != ref_file:
            rt1, rt2 = dtw(msfiles[id_ref].original_time,
                           msfiles[id_ref].bin_spectrum,
                           msfiles[_i].original_time,
                           msfiles[_i].bin_spectrum)
            loess_func_2t1 = loess_dtw(rt1, rt2)

            # map file
            new_exp = []
            for j, scan in enumerate(msfiles[_i].exp):
                old_time = msfiles[_i].original_time[j]
                new_time = loess_func_2t1(old_time).item()
                msfiles[_i].corrected_time[j] = new_time
                new_exp.append(ms1scan(scan.mz, scan.i))
            msfiles[_i].exp = new_exp
            plt.plot(msfiles[_i].corrected_time, msfiles[_i].original_time - msfiles[_i].corrected_time, color=colors[_i])

        # save file
        file_name = os.path.basename(file)
        new_path = os.path.join(new_dir, file_name)
        msfiles[_i].save_file(new_path)

    plt.xlabel('Corrected retention time (s)')
    plt.ylabel('Deviation of retention time (s)')
    plt.tight_layout()
    plt.savefig(figure_file, dpi=300)


# perform dtw and anchor alignment on a given file list (generate delta rt plot and aligned files)
def alignment_dtw_anchor(input_files, new_dir, figure_file):
    msfiles = []
    for file in input_files:
        msfiles.append(MsFile(file))

    id_ref = 0
    ref_file = input_files[id_ref]
    fig, ax = plt.subplots(figsize=(3, 2.5))
    colors = plt.cm.rainbow(np.linspace(0, 1, len(input_files)))

    roas = roa_construction(msfiles[id_ref], 0.001, 0.005, 30)
    mz_list = [i.mzmean for i in roas]
    for _i, file in enumerate(input_files):
        if file != ref_file:
            rt1, rt2 = dtw(msfiles[id_ref].original_time,
                           msfiles[id_ref].bin_spectrum,
                           msfiles[_i].original_time,
                           msfiles[_i].bin_spectrum)
            loess_func_1t2 = loess_dtw(rt2, rt1)

            eic_array = get_eic_array(msfiles[_i], mz_list, 0.01)
            f_rt, _ = local_match(roas, eic_array, loess_func_1t2)
            # map file
            new_exp = []
            for j, scan in enumerate(msfiles[_i].exp):
                old_time = msfiles[_i].original_time[j]
                new_time = f_rt(old_time).item()
                msfiles[_i].corrected_time[j] = new_time
                new_exp.append(ms1scan(scan.mz, scan.i))
            msfiles[_i].exp = new_exp
            plt.plot(msfiles[_i].corrected_time, msfiles[_i].original_time - msfiles[_i].corrected_time,
                     color=colors[_i])

        # save file
        file_name = os.path.basename(file)
        new_path = os.path.join(new_dir, file_name)
        msfiles[_i].save_file(new_path)

    plt.xlabel('Corrected retention time (s)')
    plt.ylabel('Deviation of retention time (s)')
    plt.tight_layout()
    plt.savefig(figure_file, dpi=300)


def alignment_combined(input_files, new_dir, figure_file):
    def plot_alignment(msfiles, ref_file, colors, line_style):
        for _i, file in enumerate(input_files):
            if file != ref_file:
                rt1, rt2 = dtw(msfiles[id_ref].original_time,
                               msfiles[id_ref].bin_spectrum,
                               msfiles[_i].original_time,
                               msfiles[_i].bin_spectrum)
                loess_func_2t1 = loess_dtw(rt1, rt2)

                # map file
                new_exp = []
                for j, scan in enumerate(msfiles[_i].exp):
                    old_time = msfiles[_i].original_time[j]
                    new_time = loess_func_2t1(old_time).item()
                    msfiles[_i].corrected_time[j] = new_time
                    new_exp.append(ms1scan(scan.mz, scan.i))
                msfiles[_i].exp = new_exp
                plt.plot(msfiles[_i].corrected_time, msfiles[_i].original_time - msfiles[_i].corrected_time,
                         linestyle=line_style, color=colors[_i], linewidth=0.7)

            # save file
            # file_name = os.path.basename(file)
            # new_path = os.path.join(new_dir, file_name)
            # msfiles[_i].save_file(new_path)

    msfiles = [MsFile(file) for file in input_files]
    id_ref = 0
    ref_file = input_files[id_ref]
    fig, ax = plt.subplots(figsize=(3, 2.5))
    colors = plt.cm.rainbow(np.linspace(0, 1, len(input_files)))

    # Plot alignment_only_dtw
    plot_alignment(msfiles, ref_file, colors, '-')

    # Reinitialize msfiles for the next alignment method
    msfiles = [MsFile(file) for file in input_files]
    roas = roa_construction(msfiles[id_ref], 0.001, 0.005, 30)
    mz_list = [i.mzmean for i in roas]

    def plot_alignment_anchor(msfiles, ref_file, colors, line_style):
        for _i, file in enumerate(input_files):
            if file != ref_file:
                rt1, rt2 = dtw(msfiles[id_ref].original_time,
                               msfiles[id_ref].bin_spectrum,
                               msfiles[_i].original_time,
                               msfiles[_i].bin_spectrum)
                loess_func_1t2 = loess_dtw(rt2, rt1)

                eic_array = get_eic_array(msfiles[_i], mz_list, 0.01)
                f_rt, _ = local_match(roas, eic_array, loess_func_1t2)
                # map file
                new_exp = []
                for j, scan in enumerate(msfiles[_i].exp):
                    old_time = msfiles[_i].original_time[j]
                    new_time = f_rt(old_time).item()
                    msfiles[_i].corrected_time[j] = new_time
                    new_exp.append(ms1scan(scan.mz, scan.i))
                msfiles[_i].exp = new_exp
                plt.plot(msfiles[_i].corrected_time, msfiles[_i].original_time - msfiles[_i].corrected_time,
                         linestyle=line_style, color=colors[_i], linewidth=0.7)

            # save file
            # file_name = os.path.basename(file)
            # new_path = os.path.join(new_dir, file_name)
            # msfiles[_i].save_file(new_path)

    # Plot alignment_dtw_anchor
    plot_alignment_anchor(msfiles, ref_file, colors, '--')

    plt.xlabel('Corrected retention time (s)')
    plt.ylabel('Deviation of retention time (s)')
    solid_line = plt.Line2D([], [], color='black', linestyle='-', label='DTW only')
    dashed_line = plt.Line2D([], [], color='black', linestyle='--', label='DTW+Local matching')

    # Add the proxy artists to the legend
    plt.legend(handles=[solid_line, dashed_line])
    plt.tight_layout()
    plt.savefig(figure_file, dpi=300)


def plot_anchor_type_distribution(file1_roas, file2_eics, loess_function, rtime, tic, per=68.27):
    global outlier_indices, residual_cutoff
    file1_rt, file1_mz, file2_rt, file2_mz = [], [], [], []
    file1_rt_u, file1_mz_u, file2_rt_u, file2_mz_u = [], [], [], []
    pn = 5  # Candidate peak numbers is 5

    for i in range(len(file1_roas)):
        corr = signal.correlate(file2_eics['intensity'][i], file1_roas[i].i, mode='same')
        peaks, _ = signal.find_peaks(corr)
        sorted_peaks = sorted(peaks, key=lambda h: corr[h], reverse=True)
        if len(sorted_peaks) >= 2 and corr[sorted_peaks[0]] > corr[sorted_peaks[1]] * 5:
            file1_rt.append(file1_roas[i].rt[file1_roas[i].i_max_id])
            file1_mz.append(file1_roas[i].mzmean)
            file2_rt.append(file2_eics['rt'][sorted_peaks[0]])
            file2_mz.append(file2_eics['mz'][i][sorted_peaks[0]])
        elif len(sorted_peaks) == 1:
            file1_rt.append(file1_roas[i].rt[file1_roas[i].i_max_id])
            file1_mz.append(file1_roas[i].mzmean)
            file2_rt.append(file2_eics['rt'][sorted_peaks[0]])
            file2_mz.append(file2_eics['mz'][i][sorted_peaks[0]])
        elif len(sorted_peaks) >= 2:
            file1_rt_u.append(file1_roas[i].rt[file1_roas[i].i_max_id])
            file1_mz_u.append(file1_roas[i].mzmean)
            file2_rt_u.append(
                [file2_eics['rt'][sorted_peaks[j]] if j < len(sorted_peaks) else -1000 for j in range(pn)])
            file2_mz_u.append(
                [file2_eics['mz'][i][sorted_peaks[j]] if j < len(sorted_peaks) else -1000 for j in range(pn)])
    file1_rt, file2_rt, file1_mz, file2_mz = np.array(file1_rt), np.array(file2_rt), np.array(file1_mz), np.array(
        file2_mz)
    file1_rt_u, file2_rt_u, file1_mz_u, file2_mz_u = np.array(file1_rt_u), np.array(file2_rt_u), np.array(
        file1_mz_u), np.array(file2_mz_u)

    loess_residuals = np.abs(loess_function(file1_rt) - file2_rt)
    q = 0.01
    n = len(file1_rt)
    alpha = q * (n - np.arange(n)) / n
    p68 = np.percentile(loess_residuals, per, interpolation='linear')
    rsdr = p68 * n / (n - 2)
    rank_indices = loess_residuals.argsort()
    sorted_residuals = loess_residuals[rank_indices]
    # for i in range(int(0.7 * n) - 1, n):
    for i in range(int(0.8 * n) - 1, n):
        t_value = sorted_residuals[i] / rsdr
        p_value = 2 * (1 - t.cdf(t_value, n - 2))
        if p_value < alpha[i]:
            outlier_indices = rank_indices[i:]
            residual_cutoff = sorted_residuals[i - 1]
            break
        elif i == n - 1:
            outlier_indices = []
            residual_cutoff = sorted_residuals[i]
    mask = np.isin(np.arange(n), outlier_indices, invert=True)
    rt1 = file1_rt[mask]
    rt2 = file2_rt[mask]
    mz1 = file1_mz[mask]
    mz2 = file2_mz[mask]
    if file1_rt_u.size != 0:
        pred_file1_rt_u = loess_function(file1_rt_u)
        act_minus_pred = np.abs(file2_rt_u - pred_file1_rt_u.reshape(-1, 1))
        min_id_amp = np.argmin(act_minus_pred, axis=1)
        min_value_test = np.amin(act_minus_pred, axis=1) <= residual_cutoff
        a_rt1 = file1_rt_u[np.where(min_value_test)[0]]
        a_rt2 = file2_rt_u[np.where(min_value_test)[0], min_id_amp[np.where(min_value_test)[0]]]
        a_mz1 = file1_mz_u[np.where(min_value_test)[0]]
        a_mz2 = file2_mz_u[np.where(min_value_test)[0], min_id_amp[np.where(min_value_test)[0]]]

    res_rt1, res_rt2 = np.append(rt1, a_rt1), np.append(rt2, a_rt2)
    res_mz1, res_mz2 = np.append(mz1, a_mz1), np.append(mz2, a_mz2)
    fig, ax = plt.subplots(figsize=(3, 2.5))

    ax.scatter(rt1, rt2, s=2.5, color='#e41a1c', alpha=0.8, label='Remaining type 1 anchors')
    ax.scatter(a_rt1, a_rt2, s=2.5, color='#377eb8', alpha=0.8, label='Remaining type 2 anchors')

    ax.scatter(file1_rt[~mask], file2_rt[~mask], s=5, marker='x', alpha=0.8, color='#984ea3', label='Outliers')

    t1_test = np.linspace(0, 720, 500)
    t2_pred = loess_function(t1_test)
    plt.plot(t1_test, t2_pred, color='#ff7f00', label='Profile mapping curve', linewidth=0.8)
    # plt.plot(t1_test, t2_pred - residual_cutoff, linestyle='dashed', color='red', linewidth=1)
    # plt.plot(t1_test, t2_pred + residual_cutoff, linestyle='dashed', color='red', linewidth=1)
    plt.plot(rtime, tic, color='gray', linewidth=0.8)

    plt.xlabel('Retention time (Reference file)')
    plt.ylabel('Retention time (Target file)')
    # plt.xlim((150, 300))
    # plt.ylim((150, 300))
    plt.legend()
    plt.tight_layout()
    plt.savefig('xgym_pos_anchor_type_dectection.tiff', dpi=300)


if __name__ == '__main__':
    import os
    input_dir = 'F:\\xgym_pos_500\\mzml'
    input_files_list = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.mzML')]
    input_files_list = [os.path.normpath(i) for i in input_files_list][100:110]
    input_files_list = ['F:\\xgym_pos_500\\mzml\\QC-52.mzML'] + input_files_list
    # input_files_list = ['F:\\xgym_pos_500\\mzml\\QC-52.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-98-1.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-274-4.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-229-5.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-88-4.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-177-4.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-31-2.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-104-3.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-175-2.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-49-1.mzML',
    #                     'F:\\xgym_pos_500\\mzml\\LC-297-3.mzML',]

    # alignment_only_dtw(input_files_list,
    #                    'F:\\xgym_pos_500\\dtw-new-mzml',
    #                    './alignment-dtw-only.tiff')
    # alignment_dtw_anchor(input_files_list,
    #                      'F:\\xgym_pos_500\\dtw-anchor-new-mzml',
    #                      './alignment-dtw-anchor.tiff')
    # alignment_combined(input_files_list,
    #                    'F:\\xgym_pos_500\\dtw-anchor-new-mzml',
    #                    './alignment-combined.tiff')
    sample_data = [MsFile(i) for i in input_files_list]
    rt1, rt2 = dtw(sample_data[0].original_time,
                   sample_data[0].bin_spectrum,
                   sample_data[1].original_time,
                   sample_data[1].bin_spectrum
                   )

    loess_func_2t1 = loess_dtw(rt1, rt2)
    loess_func_1t2 = loess_dtw(rt2, rt1)

    file1_roas = roa_construction(sample_data[0], int_cutoff=0.005)
    print('roa number: ', len(file1_roas))

    mzlist = [roa.mzmean for roa in file1_roas]
    file2_eic_array = get_eic_array(sample_data[1], mzlist)

    tic = []
    for scan in sample_data[0].exp:
        tic.append(scan.i.sum() / 6e6 - 30)
    plot_anchor_type_distribution(file1_roas, file2_eic_array, loess_func_1t2,
                                  sample_data[0].original_time, tic)


