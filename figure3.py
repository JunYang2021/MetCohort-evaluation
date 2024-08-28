# feature detection on ROI-matrix
# fig 3, extended fig 2, extended fig 3, fig s8

import sys

sys.path.append('../src')
from src import *
import pickle
import matplotlib.pyplot as plt
import numpy as np
from src.batch_peakdetection import row_norm_uv, gaussian_kernel, sobel_filter, non_ext_suppression, gradient_norm, peak_mask
from src.targeted_extraction import non_maxs_1d, non_mins_1d
from scipy import signal
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.ticker as ticker

font = {'family': 'Arial', 'size': 7}
plt.rc('font', **font)


def get_batch_roi(pkd_file, id):
    with open(pkd_file, 'rb') as f:
        peak_res = pickle.load(f)
        file_list = pickle.load(f)
    batch_roi = peak_res[id]['group'][0]  # BatchROI
    return batch_roi, file_list


def plot_roi_matrix(batch_roi):
    plt.imshow(batch_roi.i_array)
    plt.show()


def vis_matrix(d2m, ret_array, output_path, cmap, colorbar=None):
    # fig, ax = plt.subplots(figsize=(1.9, 2.7))
    fig, ax = plt.subplots(figsize=(1.9, 1.3))
    im = ax.imshow(d2m, aspect='auto', cmap=cmap)
    ax.set_ylabel('Sample')
    ax.set_xlabel('Retention time (s)')
    num_ticks = min(len(ret_array), 6)  # Limit the number of ticks to 10 for better readability
    tick_positions = np.linspace(0, len(ret_array) - 1, num_ticks, dtype=int)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])
    if colorbar:
        fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.15)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    # plt.show()


def vis_discrete_matrix(d2m, ret_array, output_path, colorbar=None):
    # Define a custom colormap
    cmap = ListedColormap(['#440154', '#21918c', '#fde725'])  # virdis
    bounds = [-1.5, -0.5, 0.5, 1.5]
    norm = BoundaryNorm(bounds, cmap.N)

    # fig, ax = plt.subplots(figsize=(1.9, 2.7))        # have colorbar
    fig, ax = plt.subplots(figsize=(1.9, 2.4))
    im = ax.imshow(d2m, aspect='auto', cmap=cmap, norm=norm)
    ax.set_ylabel('Sample')
    ax.set_xlabel('Retention time (s)')
    num_ticks = min(len(ret_array), 6)  # Limit the number of ticks to 10 for better readability
    tick_positions = np.linspace(0, len(ret_array) - 1, num_ticks, dtype=int)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])

    if colorbar:
        cbar = fig.colorbar(im, ax=ax, ticks=[-1, 0, 1], orientation='horizontal', pad=0.15)
        cbar.ax.set_xticklabels(['-1', '0', '1'])  # Manually set the tick labels

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)


def matrix_processing(d2m):
    d2m_1 = d2m.copy()  # original matrix
    d2m_2 = row_norm_uv(d2m_1)
    d2m_3 = signal.convolve2d(d2m_2, gaussian_kernel(7, 2), boundary='symm', mode='same')
    d2m_4 = sobel_filter(d2m_3)
    d2m_5 = non_ext_suppression(d2m_4)
    d2m_6 = gradient_norm(d2m_5)
    return d2m_1, d2m_2, d2m_3, d2m_4, d2m_5, d2m_6


def vis_matrix_processing(d2m, ret_array, output_path):
    num_ticks = min(len(ret_array), 6)  # Limit the number of ticks to 10 for better readability
    tick_positions = np.linspace(0, len(ret_array) - 1, num_ticks, dtype=int)
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(7, 6))
    ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()

    d2m_1, d2m_2, d2m_3, d2m_4, d2m_5, d2m_6 = matrix_processing(d2m)

    ax1.imshow(d2m_1, aspect='auto')
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])
    ax1.set_title('1. Original ROI-matrix')
    ax1.set_ylabel('Sample')
    ax1.set_xlabel('Retention time (s)')
    ax2.imshow(d2m_2, aspect='auto')
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])
    ax2.set_title('2. Row normalization')
    ax2.set_ylabel('Sample')
    ax2.set_xlabel('Retention time (s)')
    ax3.imshow(d2m_3, aspect='auto')
    ax3.set_xticks(tick_positions)
    ax3.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])
    ax3.set_title('3. Intensity smoothing')
    ax3.set_ylabel('Sample')
    ax3.set_xlabel('Retention time (s)')
    ax4.imshow(d2m_4, aspect='auto', cmap='gray')
    ax4.set_xticks(tick_positions)
    ax4.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])
    ax4.set_title('4. Sobel gradient')
    ax4.set_ylabel('Sample')
    ax4.set_xlabel('Retention time (s)')
    ax5.imshow(d2m_5, aspect='auto', cmap='gray')
    ax5.set_xticks(tick_positions)
    ax5.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])
    ax5.set_title('5. Non extreme suppression')
    ax5.set_ylabel('Sample')
    ax5.set_xlabel('Retention time (s)')
    ax6.imshow(d2m_6, aspect='auto', cmap='gray')
    ax6.set_xticks(tick_positions)
    ax6.set_xticklabels([f'{ret_array[i]:.0f}' for i in tick_positions])
    ax6.set_title('6. Gradient normalization')
    ax6.set_ylabel('Sample')
    ax6.set_xlabel('Retention time (s)')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)


def feature_detection_process(norm_grad_mat):
    M, N = norm_grad_mat.shape
    # find positive line
    mD_p = np.empty((M, N), dtype=np.int8)
    mS_p = np.empty((M, N), dtype=np.float16)
    mD_p[0, :] = 0
    mS_p[0, :] = norm_grad_mat[0, :]
    # find negative line
    mD_n = np.empty((M, N), dtype=np.int8)
    mS_n = np.empty((M, N), dtype=np.float16)
    mD_n[0, :] = 0
    mS_n[0, :] = norm_grad_mat[0, :]

    for i in range(1, M):
        left_p = np.pad(mS_p[i - 1, :-1], (1, 0), mode='edge') * 0.9
        up_p = mS_p[i - 1, :]
        right_p = np.pad(mS_p[i - 1, 1:], (0, 1), mode='edge') * 0.9
        mS_p[i, :] = norm_grad_mat[i, :] + np.maximum.reduce([up_p, left_p, right_p])
        directions_p = np.argmax([up_p, left_p, right_p], axis=0)
        mD_p[i, :] = np.array([0, -1, 1])[directions_p]

        left_n = np.pad(mS_n[i - 1, :-1], (1, 0), mode='edge') * 0.9
        up_n = mS_n[i - 1, :]
        right_n = np.pad(mS_n[i - 1, 1:], (0, 1), mode='edge') * 0.9
        mS_n[i, :] = norm_grad_mat[i, :] + np.minimum.reduce([up_n, left_n, right_n])
        directions_n = np.argmin([up_n, left_n, right_n], axis=0)
        mD_n[i, :] = np.array([0, -1, 1])[directions_n]

    return mS_p, mD_p, mS_n, mD_n


def border_determine(mS_p, mD_p, mS_n, mD_n):
    M, N = mS_p.shape
    res_array_list = []
    peak_res_list = []
    res_array = np.zeros((M, N))
    start_id, end_id = None, None
    find_pos = True
    non_maxs_ll = non_maxs_1d(mS_p[-1, :])
    non_mins_ll = non_mins_1d(mS_n[-1, :])
    for col in range(N):
        if find_pos:
            # if mS_p[-1, col] > M * 0.5: # M * 0.5 is threshold of line
            if non_maxs_ll[col] > M * 0.3:
                colIdx = col
                seam = [np.clip(colIdx - 2, 0, N - 1)]
                for i in range(M - 2, -1, -1):
                    colIdx = np.clip(mD_p[i + 1, colIdx] + colIdx, 0, N - 1)
                    # colIdx = np.clip(colIdx, 0, N - 1)
                    seam.append(np.clip(colIdx - 2, 0, N - 1))
                if abs(seam[-1] - seam[0]) > 10:
                    continue
                start_id = np.array(seam[::-1])
                if end_id is None or (start_id > (end_id - 3)).all():
                    find_pos = False

        else:
            # if mS_n[-1, col] < -M * 0.5:
            if non_mins_ll[col] < -M * 0.3:
                colIdx = col
                seam = [np.clip(colIdx + 2, 0, N - 1)]
                for i in range(M - 2, -1, -1):
                    colIdx = np.clip(mD_n[i + 1, colIdx] + colIdx, 0, N - 1)
                    # colIdx = np.clip(colIdx + 2, 0, N - 1)
                    seam.append(np.clip(colIdx + 2, 0, N - 1))
                if abs(seam[-1] - seam[0]) > 10:
                    continue
                end_id = np.array(seam[::-1])
                if (start_id < end_id).all():
                    find_pos = True
                    res_array[np.arange(M), start_id] = 1
                    res_array[np.arange(M), end_id] = -1
                    res_array_list.append(res_array)
                    peak_res_list.append(np.column_stack((start_id, end_id)))
                    res_array = np.zeros((M, N))

    return res_array_list, peak_res_list


def vis_integration(batch_roi, res, output_path):  # 2 features
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(2.4, 4.1))
    sample_l, points_l = batch_roi.i_array.shape
    plot_length = min(50, sample_l)
    colors = plt.cm.get_cmap('tab10', plot_length)
    start_t, end_t = batch_roi.rt_array[0], batch_roi.rt_array[-1]
    for i in range(plot_length):
        peak_start, peak_end = res[0][i]  # plot samples randomly
        # peak_start, peak_end = self.peak_pos[k][i]
        rt_range = (batch_roi.rt_array >= start_t) & (batch_roi.rt_array <= end_t)
        color = colors(i)
        ax1.plot(batch_roi.rt_array[rt_range], batch_roi.i_array[i][rt_range], linewidth=0.7, color=color)
        ax1.fill_between(batch_roi.rt_array[peak_start:peak_end + 1], batch_roi.i_array[i][peak_start: peak_end + 1],
                         color=color, alpha=0.3)
    for i in range(plot_length):
        peak_start, peak_end = res[1][i]  # plot samples randomly
        # peak_start, peak_end = self.peak_pos[k][i]
        rt_range = (batch_roi.rt_array >= start_t) & (batch_roi.rt_array <= end_t)
        color = colors(i)
        ax2.plot(batch_roi.rt_array[rt_range], batch_roi.i_array[i][rt_range], linewidth=0.7, color=color)
        ax2.fill_between(batch_roi.rt_array[peak_start:peak_end + 1], batch_roi.i_array[i][peak_start: peak_end + 1],
                         color=color, alpha=0.3)

    ax1.set_ylabel("Intensity")
    ax2.set_ylabel('Intensity')
    ax2.set_xlabel('Retention time (s)')

    formatter = ticker.ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))

    ax1.yaxis.set_major_formatter(formatter)
    ax2.yaxis.set_major_formatter(formatter)
    ax2.set_ylim((-1e5, 2e6))
    plt.tight_layout()
    # plt.show()
    plt.savefig(output_path, dpi=300)


def plot_xic(msfile, rt, mz, out_path):
    rt_list = []
    i_list = []
    for t, s in zip(msfile.original_time, msfile.exp):
        if rt[0] < t < rt[1]:
            _m, _i = s.mz, s.i
            mask1 = (_m > mz[0]) & (_m < mz[1])
            rt_list.append(t)
            i_list.append(_i[mask1].sum())
        elif t > rt[1]:
            break
    fig, ax = plt.subplots(figsize=(2.7, 2.1))
    ax.plot(rt_list, i_list, color='#1f78b4', linewidth=0.8)
    ax.set_ylabel('Intensity')  # Rotate the label and adjust the labelpad
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_xlabel('Retention time (s)')
    formatter = ticker.ScalarFormatter()
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(formatter)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)


if __name__ == '__main__':
    sample_file = [r'F:\Test_Peakmat\xgym data\mzml data\LC-3-5.mzML',
                   r'F:\Test_Peakmat\xgym data\mzml data\LC-5-4.mzML',
                   r'F:\Test_Peakmat\xgym data\mzml data\LC-20-2.mzML']
    sample_data = load_files(sample_file)
    # plot_xic(sample_data[sample_file[0]], (370, 442), (198.1742, 198.1942),
    #          r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure3\xic_198_370.tiff')
    # plot_xic(sample_data[sample_file[1]], (370, 442), (198.1742, 198.1942),
    #          r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure3\xic_198_370_lc_5_4.tiff')
    # plot_xic(sample_data[sample_file[1]], (293, 312), (991.6593, 991.6793),
    #          r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure3\xic_991_300_lc_5_4.tiff')
    plot_xic(sample_data[sample_file[1]], (293, 312), (991.6593, 991.6793),
             r'D:\Metabolomics2022\文章\peakmat文章\peakmat article 20240513\pictures\figure3\xic_991_300_lc_20_2.tiff')

    # peak 391
    # batch_roi, file_list = get_batch_roi(
    #     r'E:\yangjun\test_peakmat_v3\feature detection example\xgym_random_500_3453.pkd', 391)
    # plot_roi_matrix(batch_roi)
    # d2m = batch_roi.i_array[-200:]
    # vis_matrix(d2m, batch_roi.rt_array, './xgym_roi_matrix_391_origin.tiff')
    # vis_matrix_processing(d2m, batch_roi.rt_array, './xgym_roi_matrix_391_processing.tiff')

    # d2m_1, d2m_2, d2m_3, d2m_4, d2m_5, d2m_6 = matrix_processing(d2m)
    # vis_matrix(d2m_6, batch_roi.rt_array, './xgym_roi_matrix_391_normalized_gradient.tiff', cmap='gray')

    # mS_p, mD_p, mS_n, mD_n = feature_detection_process(d2m_6)
    # vis_matrix(mS_p, batch_roi.rt_array, './xgym_roi_matrix_391_score_positive.tiff', cmap='viridis', colorbar=True)
    # vis_discrete_matrix(mD_p, batch_roi.rt_array, './xgym_roi_matrix_391_direction_positive.tiff', colorbar=True)
    # vis_matrix(mS_n, batch_roi.rt_array, './xgym_roi_matrix_391_score_negative.tiff', cmap='viridis', colorbar=True)
    # vis_discrete_matrix(mD_n, batch_roi.rt_array, './xgym_roi_matrix_391_direction_negative.tiff', colorbar=True)

    # res, peak_res = border_determine(mS_p, mD_p, mS_n, mD_n)
    # res_array1, res_array2 = res[0], res[1]
    # vis_discrete_matrix(res_array1, batch_roi.rt_array, './xgym_roi_matrix_391_peak1.tiff')
    # vis_discrete_matrix(res_array2, batch_roi.rt_array, './xgym_roi_matrix_391_peak2.tiff')

    # vis_integration(batch_roi, peak_res, './xgym_roi_matrix_391_integration.tiff')


    # peak 3449
    batch_roi, file_list = get_batch_roi(
        r'E:\yangjun\test_peakmat_v3\feature detection example\xgym_random_500_3453.pkd', 3449)
    # plot_roi_matrix(batch_roi)
    d2m = batch_roi.i_array[-100:]
    # vis_matrix(d2m, batch_roi.rt_array, './xgym_roi_matrix_3449_origin.tiff', None)
    # vis_matrix_processing(d2m, batch_roi.rt_array, './xgym_roi_matrix_3449_processing.tiff')

    d2m_1, d2m_2, d2m_3, d2m_4, d2m_5, d2m_6 = matrix_processing(d2m)
    # vis_matrix(d2m_6, batch_roi.rt_array, './xgym_roi_matrix_3449_normalized_gradient.tiff', cmap='gray')
    mS_p, mD_p, mS_n, mD_n = feature_detection_process(d2m_6)
    res, peak_res = border_determine(mS_p, mD_p, mS_n, mD_n)
    mask_g = np.where(peak_mask(d2m_5, peak_res, d=3), d2m_5, 0)
    new_d2m_6 = gradient_norm(mask_g)
    # vis_matrix(new_d2m_6, batch_roi.rt_array, './xgym_roi_matrix_3449_normalized_gradient_after_mask.tiff', cmap='gray')
    mS_p, mD_p, mS_n, mD_n = feature_detection_process(new_d2m_6)
    res, new_peak_res = border_determine(mS_p, mD_p, mS_n, mD_n)
    peak_res = peak_res + new_peak_res
    vis_integration(batch_roi, peak_res, './xgym_roi_matrix_3449_integration_new.tiff')