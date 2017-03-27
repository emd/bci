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

    vibration_subtracted - bool
        If True and `self.beam` is CO2, the vibrational contributions
        to the phase data `self.x` have been removed.

    Methods:
    --------
    t - returns retrieved signal time-base, array-like, (`N`,)
        The time-base is generated on the fly as needed and is not stored
        as an object property; this helps save memory and processing time,
        as we do not typically look at the raw signal vs. time.
        [t] = s

    '''
    def __init__(self, shot, chord='V2', beam='CO2',
                 tlim=[-0.05, 5.2], vibration_subtracted=False):
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

            The specified `tlim` values will always bracket the retrieved
            signal. That is, if `tlim[0]` does not correspond to an exact
            digitization time, then the initial time returned (`Signal.t0`)
            will be the closest digitization time *greater* than `tlim[0]`.
            Similarly, if `tlim[1]` does not correspond to an exact
            digitization time, then the final time (`Signal.t[-1]`) will be
            the closest digitization time *less* than `tlim[1]`. Further,
            if `tlim[0]` is less than the initial digitization time,
            the retrieved signal will begin with the initial digitized point.
            Similarly, if `tlim[1]` exceeds the final digitization time,
            the retrieved signal will end with the final digitized point.

            If `tlim` is not `None` and it is *not* a length-two array,
            a ValueError is raised.

            [tlim] = s

        vibration_subtracted - bool
            If True *and* `beam` is CO2, use vibration-subtracted phase data.

        '''
        self.shot = shot

        if str.upper(chord) in set(['V1', 'V2', 'V3', 'R0']):
            self.chord = str.upper(chord)
        else:
            raise ValueError('`chord` may be V1, V2, V3, or R0.')

        if str.upper(beam) == 'CO2':
            self.beam = 'CO2'
        elif str.upper(beam) == 'HENE':
            self.beam = 'HeNe'
        else:
            raise ValueError('`beam` may be CO2 or HeNe.')

        if vibration_subtracted and (self.beam == 'CO2'):
            self.vibration_subtracted = True
        else:
            self.vibration_subtracted = False

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
        sig = np.zeros(_Npts_per_window * len(windows))

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

        if self.vibration_subtracted and (self.beam == 'CO2'):
            phase_datatype = 'vibration-subtracted'
        else:
            phase_datatype = 'raw'

        print ('\nLoading %s %s phase data (%s)'
               % (self.chord, self.beam, phase_datatype))

        # Cycle through each window to iteratively build up signal
        for i, window in enumerate(windows):
            print 'Window %i (%i of %i)' % (window, i + 1, len(windows))

            # Determine whether to retrieve raw or vibration-subtracted phase
            if self.vibration_subtracted and (self.beam == 'CO2'):
                # [line-integrated density] = m / cm^3
                node = tree.getNode('\DEN%s_UF_%i' % (self.chord, window))
            else:
                # [phase] = radians
                node = tree.getNode(
                    '\PL%i%s_UF_%i' % (beam_number, self.chord, window))

            # Insert retrieved signal into appropriate location
            sl = slice(_Npts_per_window * i, (_Npts_per_window * (i + 1)))
            sig[sl] = node.getData().data()

        # Crop signal to desired time window
        t0, sig = _crop(sig, tlim)

        # Convert from double-pass to single-pass measurement
        sig /= 2

        # Convert vibration-subtracted, line-integrated density
        # to corresponding phase, if needed.
        if self.vibration_subtracted and (self.beam == 'CO2'):
            re = 2.818e-15      # classical electron radius, [re] = m
            lambda0 = 10.6e-6   # CO2 wavelength, [lambda0] = m
            cm_per_m = 100
            sig *= (re * lambda0 * (cm_per_m ** 3))

        return t0, sig

    def t(self):
        'Get times for points in `self.x`.'
        return self.t0 + (np.arange(len(self.x)) / self.Fs)


def _closest_digitized_point(t, mode=np.round):
    '''Get global index of digitized point closest to time `t`,
    subject to constraint provided by `mode`.

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

    mode - callable
        May be {`np.round`, `np.floor`, `np.ceil`},
        otherwise a ValueError is raised.

        - `np.round`: get closest digitized point to time `t`
        - `np.floor`: get digitized point that is closest to
          *and* less-than or equal to time `t`
        - `np.ceil`: get digitized point that is closest to
          *and* greater-than or equal to time `t`

    Returns:
    --------
    pt - int
        The global index of digitized point closest to time `t`,
        subject to constraint provided by `mode`.

    '''
    if mode in set([np.round, np.floor, np.ceil]):
        pt = np.int(mode(_Fs * (t - _trigger_time)))
    else:
        raise ValueError('`mode` may be {`np.round`, `np.floor`, `np.ceil`}')

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
    pt_min = _closest_digitized_point(tlim[0], mode=np.ceil)
    pt_max = _closest_digitized_point(tlim[1], mode=np.floor)

    wlo = pt_min // _Npts_per_window
    whi = pt_max // _Npts_per_window

    return np.arange(wlo, whi + 1)


def _crop(sig, tlim):
    'Crop `sig` to time window `tlim`.'
    # First, determine start and stop indices in global record
    gstart = _closest_digitized_point(tlim[0], mode=np.ceil)
    gstop = _closest_digitized_point(tlim[1], mode=np.floor)

    # Convert global indices to "local" indices relevant for slicing `sig`
    lstart = gstart % _Npts_per_window
    lstop = lstart + (gstop - gstart)

    # Determine time closest to `tlim[0]`
    t0 = _trigger_time + (gstart / _Fs)

    return t0, sig[lstart:(lstop + 1)]


def _plasma_induced_phase(
        ph1, ph2, lambda1=10.6e-6, lambda2=0.633e-6):
    '''Get plasma-induced phase corresponding to raw phase measurement
    `ph1` by using phase measurements at another wavelength to
    subtract the vibrational contributions from `ph1`.

    Parameters:
    -----------
    ph1 - array_like, (`N`,)
        The raw phase measurement at probe wavelength `lambda1`.
        Vibrational contributions will be removed from `ph1` such that
        the returned phase consists only of the plasma-induced
        contributions to `ph1`.
        [ph1] = radian

    ph2 - array_like, (`N`,)
        The raw phase measurement at probe wavelength `lambda2`.
        This measurement is used to remove the vibrational contributions
        from `ph1`.
        [ph2] = radian

    lambda1 - float
        Probe wavelength of the primary beam.
        [lambda1] = m

    lambda2 - float
        Probe wavelength of the secondary beam.
        [lambda2] = m

    Returns:
    --------
    ph1_plasma - array_like, (`N`,)
        The plasma-induced phase contributions to `ph1`.
        [ph1_plasma] = radian

    '''
    num = lambda1 * ((lambda1 * ph1) - (lambda2 * ph2))
    den = (lambda1 ** 2) - (lambda2 ** 2)

    return num / den
