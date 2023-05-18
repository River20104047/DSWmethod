"""
By Zijiang Yang 2023/05/19
Matlab translated Python equivalent...hopefully

There are some differences between Matlab functions and Python functions
Some parameters were adjusted:
1. sensitivity sn = 2
I am still working on improving Python version
"""


import numpy as np
from scipy.interpolate import pchip_interpolate
from scipy.stats import gaussian_kde
import pandas as pd
import matplotlib.pyplot as plt


def dsw_method(X, Y, lbws1, check=True):
    """
    This function applies various signal processing techniques to a given spectrum.

    Parameters:
    X (np.array): The independent variable (e.g., wave number).
    Y (np.array): The dependent variable (e.g., signal amplitude).
    lbws1 (int): The small window size.
    check (bool): A flag to create a plot of the processed spectrum. Defaults to True.

    Returns:
    X (np.array): Processed independent variable.
    Ybc (np.array): Processed dependent variable.
    SNR (list): A list containing values related to the Signal-to-Noise Ratio.
    """

    def compute_envelope(x, y, es, envelope_type):
        """
        Compute the upper or lower envelope of a signal.

        Parameters:
        x (np.array): The independent variable (e.g., wave number).
        y (np.array): The dependent variable (e.g., signal amplitude).
        es (int): The size of the moving window for local min/max calculation.
        envelope_type (str): Either 'min' or 'max', specifies which envelope to compute.

        Returns:
        yev (np.array): The computed envelope.
        """
        assert envelope_type in ['min', 'max'], "envelope_type must be either 'min' or 'max'"
        x_extrema_points = []
        y_extrema_points = []
        for i in range(0, len(x), es):
            window = y[i:i + es]
            if envelope_type == 'max':
                extremum_val = np.max(window)
                extremum_val_idx = i + np.argmax(window)
            else:  # envelope_type == 'min'
                extremum_val = np.min(window)
                extremum_val_idx = i + np.argmin(window)
            x_extrema_points.append(x[extremum_val_idx])
            y_extrema_points.append(extremum_val)
        # Conduct shape-preserving interpolation based on obtained extrema signals
        yev = pchip_interpolate(x_extrema_points, y_extrema_points, x)
        return yev

    # Internal adjustable parameters (I)
    lbws2 = 5 * lbws1  # Large window size
    ex = 1.3143  # Expansion factor
    fn = 1.5  # Skewness correction factor (rough)
    sn = 2  # Sensitivity of peak detection (it was 4 in Matlab)
    fs = 1  # Smoothing factor
    pr = lbws1  # Smoothing radium when combining baselines
    es = lbws1  # Window size and step size for baseline correction

    # Internal empirical parameters (II)
    fup = 3.777 * lbws1 ** -0.5296 + 0.8164
    flw = 0.9594 * lbws1 ** -0.6008 + 0.9407
    fsk = 0.2246 * lbws1 ** -0.1126 + 0.2369

    # 0. Negative signal correction
    Y = Y - min(Y)

    # 1. Make sure X is in ascending order
    if X[0] > X[-1]:
        X = np.flip(X)
        Y = np.flip(Y)

    # 2. Prepare for X, Y, and XY
    XY = np.column_stack((X, Y))

    # 3. Compute the upper and lower envelopes
    Yup = compute_envelope(X, Y, es, 'max')
    Ylw = compute_envelope(X, Y, es, 'min')

    # 4. Difference between lower bound and upper bound
    dYe = Yup - Ylw

    # 5. Remove the top 99% for faster calculation
    threshold = np.percentile(dYe, 99)
    dYe = dYe[dYe <= threshold]

    # 6. Calculate kernel density function
    density_function = gaussian_kde(dYe)
    min_dYe = np.min(dYe)
    max_dYe = np.max(dYe)
    Xi = np.linspace(min_dYe, max_dYe, 5000)
    F = density_function(Xi)

    # 7. Noise range estimation
    max_F_index = np.argmax(F)
    Rn = Xi[max_F_index]

    # 8. Pre-baseline correction
    BC = compute_envelope(X, Y, lbws1, 'min')  # estimated baseline
    BC0 = BC + 0.5 * Rn  # unbiased estimate of baseline
    Ybc0 = Y - BC0  # unbiased estimate of baseline corrected spectrum

    # 9. Pre-noise removal
    Ybn0 = Ybc0 - fn * Rn  # fn * Rn is the radius of the noise removal range
    index0 = Ybn0 > 0
    Ybn0 = Ybn0 * index0

    # 10. Pre-output
    XYbc0 = np.column_stack((X, Ybc0))  # Combine X and Ybc into XYbc matrix
    XYbn0 = np.column_stack((X, Ybn0))  # Combine X and Ybn into XYbn matrix

    # 11a. Double baseline 1: small window size for noise
    Ybc1 = compute_envelope(X, Y, lbws1, 'min')
    BC1 = Ybc1 + flw * 0.5 * Rn

    # 11b. Double baseline 2: large window size for noise
    Ybc2 = compute_envelope(X, Y, lbws2, 'min')
    BC2 = Ybc2 + flw * 0.5 * Rn

    # 12. Selection index based on candidate peaks
    index0 = XYbn0[:, 1] > 0
    window_size = 4 * int(fs * pr / abs(X[1] - X[0]))  # Calculate window size based on fs and pr　()
    weights = np.ones(window_size)
    s = pd.Series(index0)  # convert it to a pandas Series
    index1 = s.rolling(window_size, center=True).mean()  # calculate the moving average
    index2 = index1 > sn / pr
    index2 = index2 + 0

    # 13. Selected baselines
    window_size = int(lbws2 / abs(X[1] - X[0]))
    s = pd.Series(index2)  # convert it to a pandas Series
    index = s.rolling(window_size, center=True).mean()  # calculate the moving average
    BCnp = BC1 * (1 - index)  # Selected baseline at non-peak region
    BCpk1 = BC1 * index  # Small BC at peak region
    BCpk2 = BC2 * index  # Large BC at peak region
    BCpk = np.minimum(BCpk1, BCpk2)
    BC_final = BCnp + BCpk

    # 14. Baseline corrected spectrum
    Ybc = Y - BC_final
    XYbc = np.column_stack((X, Ybc))

    # 15. SNR estimation
    sigma = fsk * (Rn * (fup + flw) / (1.96 * 2))  # Here, Rn * (fup + flw) is 95% CI
    H_peak = np.max(Ybc)
    snr = H_peak / sigma
    SNR = [H_peak, sigma, snr]

    # -1. Post hoc process
    H_peak_rounded = round(H_peak, 2)
    sigma_rounded = round(sigma, 2)
    snr_rounded = round(snr, 2)
    print("Peak height:", H_peak_rounded)
    print("sigma:", sigma_rounded)
    print("SNR:", snr_rounded)

    # -2. Create a figure with subplots
    if check:
        fig, axs = plt.subplots(2, 1, figsize=(10, 10))
        index_fill = index0 * max(Y)
        axs[0].fill_between(X, 0, index_fill, where=(index_fill > 0), color='#ADD8E6', alpha=0.5)
        axs[0].plot(X, BC1, color='#FF7C80', label='Small window', linewidth=0.5)
        axs[0].plot(X, BC2, color='#FF7C80', label='Large window', linewidth=0.5)
        axs[0].plot(X, Y, color='#6699FF', label='Original spectrum', linewidth=1)
        axs[0].plot(X, BC_final, color='#FF7C80', label='Estimated baseline', linewidth=1.5)
        axs[0].grid(False)
        axs[0].set_xlim(400, 4000)
        axs[0].set_ylim(0, 1.1 * max(Y))
        axs[0].invert_xaxis()
        axs[1].plot(X, Ybc, color='#00CC99', linewidth=1)
        axs[1].set_xlabel('Wavenumber (cm-1)')
        axs[1].set_ylabel('Intensity')
        axs[1].set_title('Corrected Spectrum')
        axs[1].grid(False)
        axs[1].set_xlim(400, 4000)
        axs[1].invert_xaxis()
        plt.subplots_adjust(hspace=0.5)
        plt.show()
    return X, Ybc, SNR
