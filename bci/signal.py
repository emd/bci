import numpy as np
import MDSplus as mds


# Sampling rate
# [Fs] = samples / s
Fs = (5. / 3) * 1e6

# Nominal trigger time
# [t0] = s
t0 = -1.4986927509307861

# Parameters of "ultrafast" windows
Nwindows = 9
Npts_per_window = 2 ** 21
Npts_total = (Nwindows * Npts_per_window)

# Prefactor to convert line-integrated density to corresponding phase
# for a CO2 beam
# [re_lambda0] = m^2
re_lambda0 = (2.8e-15) * (10.6e-6)


def closest_digitized_point(t):
    '''Get digitized point that most closely corresponds to time `t`,
    where [t] = s.

    '''
    pt = np.int(np.round(Fs * (t - t0)))

    # Ensure point does not exceed bounds of digitizer record
    if pt < 0:
        pt = 0
    elif pt >= Npts_total:
        pt = Npts_total - 1

    return pt


def windows(tlim):
    '''Get list of BCI windows corresponding to `tlim`.

    Parameters:
    -----------
    tlim - list, of length 2
        Lower and upper bounds of desired record in seconds.
        [tlim] = s

    Returns:
    --------
    w - array_like
        List of BCI windows corresponding to `tlim`.

    '''
    # Ensure `tlim` is ordered chronologically
    tlim = np.sort(tlim)

    # Find global position (i.e. if all windows were concatenated)
    # of the lower and upper temporal bounds
    pt_min = closest_digitized_point(tlim[0])
    pt_max = closest_digitized_point(tlim[1])

    wlo = pt_min // Npts_per_window
    whi = pt_max // Npts_per_window

    return np.arange(wlo, whi + 1)


def phase_co2(shot, tlim, chord='V2', from_density=False):
    '''Returns BCI CO2 phase along `chord` for `tlim` in `shot`.

    Parameters:
    -----------
    shot - int
        DIII-D shot number

    tlim - list, of length 2
        Requested lower and upper bounds of retrieved phase.
        [tlim] = s

    chord - string
        BCI chord: {'V1', 'V2', 'V3', 'R0'}

    from_density - bool
        Convert vibration-subtracted line-integrated density
        to corresponding phase if True (this can help remove
        low-frequency vibrational contributions in correlation
        with 285 interferometer).

    Returns:
    --------
    (tstart, ph) - tuple, where

    tstart - int
        Time (in seconds) corresponding to first point of returned phase,
        i.e. the time corresponding to `ph[0]`
        [tstart] = s

    ph - array_like
        The retrieved CO2 phase.
        [ph] = rad

    '''
    # Ensure `tlim` is ordered chronologically
    tlim = np.sort(tlim)

    ws = windows(tlim)
    tree = mds.Tree('bci', shot, 'ReadOnly')

    print '\nLoading %s CO2 phase data' % chord

    # Initialize array to store (potentially) concatenated phase data
    ph = np.zeros(Npts_per_window * len(ws))

    for i, w in enumerate(ws):
        print 'Window %i of %i' % (i + 1, len(ws))

        if from_density:
            # [line-integrated density] = m / cm^3
            node = tree.getNode('\den%s_uf_%i' % (chord, w))
        else:
            node = tree.getNode('\pl1%s_uf_%i' % (chord, w))

        sl = slice(Npts_per_window * i, (Npts_per_window * (i + 1)))
        ph[sl] = node.getData().data()

        # Convert vibration-subtracted, line-integrated density
        # to corresponding phase, if needed.
        # (The (100 ** 3) is to convert from cm^{-3} to m^{-3}).
        if from_density:
            ph[sl] *= re_lambda0 * (100 ** 3)

    # Crop returned phase to desired time window

    # First, determine position of start and stop in global record
    pt_start = closest_digitized_point(tlim[0])
    pt_stop = closest_digitized_point(tlim[1])

    # Convert positions in global record to indices of `ph`
    ind_start = pt_start % Npts_per_window
    ind_stop = ((len(ws) - 1) * Npts_per_window) + (pt_stop % Npts_per_window)

    # ... then crop
    ph = ph[ind_start:(ind_stop + 1)]

    # Determine time corresponding to `ph[0]`
    tstart = t0 + (closest_digitized_point(tlim[0]) / Fs)

    return tstart, ph


def timebase(t0, Npts):
    return t0 + (np.arange(Npts) / Fs)
