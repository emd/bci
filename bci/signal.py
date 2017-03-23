'''This module implements a class for retrieving signals from
the DIII-D bi-color interferometer (BCI) system.

'''


import numpy as np
import MDSplus as mds


# Nominal trigger time
# [_trigger_time] = s
_trigger_time = -1.4986927509307861

# Nominal sample rate
# [_Fs] = samples / s
_Fs = (5. / 3) * 1e6

# Parameters of "ultrafast" windows
_Nwindows = 9
_Npts_per_window = 2 ** 21
_Npts_total = (_Nwindows * _Npts_per_window)

# Prefactor to convert line-integrated density to corresponding phase
# for a CO2 beam
# [re_lambda0] = m^2
# re_lambda0 = (2.8e-15) * (10.6e-6)


class Signal(object):
    '''An object corresponding to the signal retrieved from the BCI system.

    Attributes:
    -----------
    shot - int
        DIII-D shot number.

    chord - string
        The interferometer chord.

    beam - string
        The type of probe beam.

    x - array-like, (`N`,)
        The retrieved, single-pass phase signal. While all of the BCI chords
        employ a Michelson configuration (and thus make a double pass through
        the plasma), this class compensates for the double-pass measurement
        for easier comparison to the MIT PCI system, which employs a
        Mach-Zehnder (i.e. single-pass) configuration.
        [x] = radian

    Fs - float
        The sampling rate.
        [Fs] = samples / second

    t0 - float
        The time corresponding to the first retrieved point in the signal;
        that is, if x(t) corresponds to the continuous signal being sampled,
        then `Signal.x[0]` = x(t0)
        [t0] = s

    Methods:
    --------
    t - returns retrieved signal time-base, array-like, (`N`,)
        The time-base is generated on the fly as needed and is not stored
        as an object property; this helps save memory and processing time,
        as we do not typically look at the raw signal vs. time.
        [t] = s

    '''
    def __init__(self, shot, chord='V2', beam='CO2', tlim=[-0.05, 5.2]):
        '''Create an instance of the `Signal` class.

        Input parameters:
        -----------------
        shot - int
            DIII-D shot number.

        chord - string
            The interferometer chord. Valid values are 'V1', 'V2',
            'V3', and 'R0'; a ValueError is raised for other values.

        beam - string
            The type of probe beam. Valid values are 'CO2' and 'HeNe';
            a ValueError is raised for other values.

        tlim - array_like, (2,)
            The lower and upper limits in time for which the signal
            will be retrieved. The default values correspond to
            the default digitization limits of the mitpci system.

            If `tlim` is not `None` and it is *not* a length-two array,
            a ValueError is raised.

            [tlim] = s

        '''
        self.shot = shot

        if str.upper(chord) in set(['V1', 'V2', 'V3', 'R0']):
            self.chord = chord
        else:
            raise ValueError('`chord` may be V1, V2, V3, or R0.')

        if str.upper(beam) == 'CO2':
            self.beam = 'CO2'
        elif str.upper(beam) == 'HENE':
            self.beam = 'HeNe'
        else:
            raise ValueError('`beam` may be CO2 or HeNe.')

        # Sampling rate
        self.Fs = _Fs

        # Phase signal between `tlim`
        self.t0, self.x = self._getSignal(tlim)

    def _getSignal(self, tlim):
        'For window `tlim`, get initial time and (phase) signal.'
        if tlim is not None:
            # Ensure limits in time are correctly sized and sorted
            if len(tlim) != 2:
                raise ValueError('`tlim` must be an iterable of length 2.')
            else:
                tlim = np.sort(tlim)

        # The BCI record for each chord and beam in any given shot
        # is split across `_Nwindows` windows; determine which windows
        # to load data from.
        windows = _windows(tlim)

        # Initialize array to store (potentially) concatenated phase data
        ph = np.zeros(_Npts_per_window * len(windows))

        # Open MDSplus tree
        tree = mds.Tree('bci', self.shot, 'ReadOnly')

        # The MDSplus node for each beam type is specified
        # by a number rather than a string
        if self.beam == 'CO2':
            beam_number = 1
        elif self.beam == 'HeNe':
            beam_number = 2
        else:
            raise ValueError('%s is not a valid beam type' % self.beam)

        print '\nLoading %s %s phase data' % (self.chord, self.beam)

        # Cycle through each window to iteratively build up signal
        for i, window in enumerate(windows):
            print 'Window %i (%i of %i)' % (window, i + 1, len(windows))

            #  if from_density:
            #      # [line-integrated density] = m / cm^3
            #      node = tree.getNode('\den%s_uf_%i' % (chord, w))
            #  else:
            node = tree.getNode(
                '\PL%i%s_UF_%i' % (beam_number, self.chord, window))

            sl = slice(_Npts_per_window * i, (_Npts_per_window * (i + 1)))
            ph[sl] = node.getData().data()

            # # Convert vibration-subtracted, line-integrated density
            # # to corresponding phase, if needed.
            # # (The (100 ** 3) is to convert from cm^{-3} to m^{-3}).
            # if from_density:
            #     ph[sl] *= re_lambda0 * (100 ** 3)

        # Crop signal to desired time window
        t0, ph = _crop(ph, tlim)

        # Convert from double-pass to single-pass measurement
        ph /= 2

        return t0, ph

    def t(self):
        'Get times for points in `self.x`.'
        return self.t0 + (np.arange(len(self.x)) / self.Fs)


def _closest_digitized_point(t):
    '''Get global index of digitized point closest to time `t`.

    Here, "global index" refers to the fact that the BCI record
    for each chord and beam in any given shot is split across
    `_Nwindows` windows. The global index specifies the position
    of a given point in the single array that results from
    sequentially concatenating the data in all of the `_Nwindows`
    windows.

    Parameters:
    -----------
    t - float
        Time.
        [t] = s

    Returns:
    --------
    pt - int
        The global index of digitized point closest to time `t`.

    '''
    pt = np.int(np.round(_Fs * (t - _trigger_time)))

    # Ensure point does not exceed bounds of digitizer record
    if pt < 0:
        pt = 0
    elif pt >= _Npts_total:
        pt = _Npts_total - 1

    return pt


def _windows(tlim):
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
    # Find global position (i.e. if all windows were concatenated)
    # of the lower and upper temporal bounds
    pt_min = _closest_digitized_point(tlim[0])
    pt_max = _closest_digitized_point(tlim[1])

    wlo = pt_min // _Npts_per_window
    whi = pt_max // _Npts_per_window

    return np.arange(wlo, whi + 1)


def _crop(sig, tlim):
    'Crop `sig` to time window `tlim`.'
    # First, determine start and stop indices in global record
    gstart = _closest_digitized_point(tlim[0])
    gstop = _closest_digitized_point(tlim[1])

    # Convert global indices to "local" indices relevant for slicing `sig`
    lstart = gstart % _Npts_per_window
    lstop = lstart + (gstop - gstart)

    # Determine time closest to `tlim[0]`
    t0 = _trigger_time + (gstart / _Fs)

    return t0, sig[lstart:(lstop + 1)]
# 
# 
# def phase_co2(shot, tlim, chord='V2', from_density=False):
#     '''Returns BCI CO2 phase along `chord` for `tlim` in `shot`.
# 
#     Parameters:
#     -----------
#     shot - int
#         DIII-D shot number
# 
#     tlim - list, of length 2
#         Requested lower and upper bounds of retrieved phase.
#         [tlim] = s
# 
#     chord - string
#         BCI chord: {'V1', 'V2', 'V3', 'R0'}
# 
#     from_density - bool
#         Convert vibration-subtracted line-integrated density
#         to corresponding phase if True (this can help remove
#         low-frequency vibrational contributions in correlation
#         with 285 interferometer).
# 
#     Returns:
#     --------
#     (tstart, ph) - tuple, where
# 
#     tstart - int
#         Time (in seconds) corresponding to first point of returned phase,
#         i.e. the time corresponding to `ph[0]`
#         [tstart] = s
# 
#     ph - array_like
#         The retrieved CO2 phase.
#         [ph] = rad
# 
#     '''
#     # Ensure `tlim` is ordered chronologically
#     tlim = np.sort(tlim)
# 
#     ws = windows(tlim)
#     tree = mds.Tree('bci', shot, 'ReadOnly')
# 
#     print '\nLoading %s CO2 phase data' % chord
# 
#     # Initialize array to store (potentially) concatenated phase data
#     ph = np.zeros(_Npts_per_window * len(ws))
# 
#     for i, w in enumerate(ws):
#         print 'Window %i of %i' % (i + 1, len(ws))
# 
#         if from_density:
#             # [line-integrated density] = m / cm^3
#             node = tree.getNode('\den%s_uf_%i' % (chord, w))
#         else:
#             node = tree.getNode('\pl1%s_uf_%i' % (chord, w))
# 
#         sl = slice(_Npts_per_window * i, (_Npts_per_window * (i + 1)))
#         ph[sl] = node.getData().data()
# 
#         # Convert vibration-subtracted, line-integrated density
#         # to corresponding phase, if needed.
#         # (The (100 ** 3) is to convert from cm^{-3} to m^{-3}).
#         if from_density:
#             ph[sl] *= re_lambda0 * (100 ** 3)
# 
#     # Crop returned phase to desired time window
# 
#     # First, determine position of start and stop in global record
#     pt_start = closest_digitized_point(tlim[0])
#     pt_stop = closest_digitized_point(tlim[1])
# 
#     # Convert positions in global record to indices of `ph`
#     ind_start = pt_start % _Npts_per_window
#     ind_stop = ((len(ws) - 1) * _Npts_per_window) + (pt_stop % _Npts_per_window)
# 
#     # ... then crop
#     ph = ph[ind_start:(ind_stop + 1)]
# 
#     # Determine time corresponding to `ph[0]`
#     tstart = t0 + (closest_digitized_point(tlim[0]) / _Fs)
# 
#     return tstart, ph
