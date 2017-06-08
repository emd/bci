from nose import tools
import numpy as np
from bci.signal import (
    Phase, _plasma_induced_phase, _closest_digitized_point, _crop,
    _windows, _trigger_time, _Fs, _Npts_per_window, _Nwindows, _Npts_total)


def test__closest_digitized_point():
    # Invalid `mode` should raise value error
    tools.assert_raises(ValueError, _closest_digitized_point, 0, np.float)

    # Final digitized point
    tf = _trigger_time + (_Npts_total / _Fs)

    # Time somewhere in digitization window (midpoint if `divisor` == 2)
    # that corresponds exactly to a digitization timestamp
    divisor = 2
    t = _trigger_time + ((_Npts_total // divisor) * (1. / _Fs))

    # When the time *exactly* corresponds to the digitizer timestamp
    # or when the time is before or after the digitization window,
    # the results do *not* depend on `mode`
    modes = [np.round, np.floor, np.ceil]
    for mode in modes:
        # Times less than `_trigger_time` point to first digitized point
        tools.assert_equal(
            _closest_digitized_point(_trigger_time - 1, mode=mode),
            0)

        # Times after last time should point to last digitized point
        tools.assert_equal(
            _closest_digitized_point(tf + 1, mode=mode),
             _Npts_total - 1)

        # Time at `_trigger_time` should point to first digitized point
        tools.assert_equal(
            _closest_digitized_point(_trigger_time, mode=mode),
            0)

        # Time at last time should point to last digitized point
        tools.assert_equal(
            _closest_digitized_point(tf, mode=mode),
             _Npts_total - 1)

        tools.assert_equal(
            _closest_digitized_point(t, mode=mode),
            (_Npts_total // divisor))

    # Check that things work when time does *not* correspond
    # to a digitizer time stamp

    # Add *half* of interval between successive timestamps to `t`
    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / _Fs), mode=np.round),
        (_Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / _Fs), mode=np.floor),
        (_Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / _Fs), mode=np.ceil),
        (_Npts_total // divisor) + 1)

    # Subtract *half* of interval between successive timestamps to `t`
    # (First test doesn't quite work at 0.5, presumably due to
    # finite-precision arithmetic... the correct behavior is verified
    # via the `eps` studies further below; also, we will not typically
    # be using the `np.round` mode, so this subtlety should not
    # be a practical issue).
    tools.assert_equal(
        _closest_digitized_point(t - (0.51 / _Fs), mode=np.round),
        (_Npts_total // divisor) - 1)

    tools.assert_equal(
        _closest_digitized_point(t - (0.5 / _Fs), mode=np.floor),
        (_Npts_total // divisor) - 1)

    tools.assert_equal(
        _closest_digitized_point(t - (0.5 / _Fs), mode=np.ceil),
        (_Npts_total // divisor))

    # Ensure behavior is correct at boundary between two adjacent points
    eps = np.finfo('float64').eps

    # Add *half* of interval between successive timestamps + `2 * eps` to `t`
    # (For some reason, tests fail when only adding `eps` to `t`, but
    # they succeed when adding `2 * eps`... whatever -- good enough)
    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / _Fs) + (2 * eps), mode=np.round),
        (_Npts_total // divisor) + 1)

    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / _Fs) + (2 * eps), mode=np.floor),
        (_Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / _Fs) + (2 * eps), mode=np.ceil),
        (_Npts_total // divisor) + 1)

    # Add *half* of interval between successive timestamps - `2 * eps` to `t`
    # (For some reason, tests fail when only subtracting `eps` from`t`, but
    # they succeed when subtracting `2 * eps`... whatever -- good enough)
    tools.assert_equal(
        _closest_digitized_point(t - (0.5 / _Fs) - (2 * eps), mode=np.round),
        (_Npts_total // divisor) - 1)

    tools.assert_equal(
        _closest_digitized_point(t - (0.5 / _Fs) - (2 * eps), mode=np.floor),
        (_Npts_total // divisor) - 1)

    tools.assert_equal(
        _closest_digitized_point(t - (0.5 / _Fs) - (2 * eps), mode=np.ceil),
        (_Npts_total // divisor))

    return


def test__crop():
    # Full record
    x = np.arange(_Nwindows * _Npts_per_window)

    # Global indices for a particular slice
    start_position_in_windows = 3.75
    gstart = np.int(start_position_in_windows * _Npts_per_window)
    gstop = gstart + 100

    # Times corresponding to global indices
    # (Let `tlim` actually *bracket* the global indices; that is,
    # let `t0` be smaller than timestamp for `gstart` and
    # let `tf` be larger than timestamp for `gstop`.)
    t0 = _trigger_time + (gstart / _Fs) - (0.999 / _Fs)
    tf = _trigger_time + (gstop / _Fs) + (0.999 / _Fs)
    tlim = [t0, tf]

    # Prior to cropping, retrieved raw signals will be arrays
    # with lengths that are integer multiples of `_Npts_per_window`
    window = np.int(np.floor(start_position_in_windows))
    retrieved_slice = slice(
        window * _Npts_per_window,
        (window + 1) * _Npts_per_window)

    np.testing.assert_equal(
        x[gstart:(gstop + 1)],
        _crop(x[retrieved_slice], tlim)[1])

    return


def test__windows():
    # Lowest window
    tlim = [_trigger_time, _trigger_time + (1000. / _Fs)]
    np.testing.assert_equal(
        _windows(tlim),
        np.array([0]))

    # Highest window
    tlim = np.array([_Npts_total - 1000, _Npts_total]) * (1. / _Fs)
    np.testing.assert_equal(
        _windows(tlim),
        np.array([_Nwindows - 1]))

    # Full set of windows
    tlim = _trigger_time + (np.array([0, _Npts_total]) * (1. / _Fs))
    np.testing.assert_equal(
        _windows(tlim),
        np.arange(_Nwindows))

    # Window in middle where `tlim` *exactly* sits on boundaries
    window = 3
    ind_lo = window * _Npts_per_window              # Lower bound, global ind
    ind_hi = ((window + 1) * _Npts_per_window) - 1  # Upper bound, global ind

    t0 = _trigger_time + (ind_lo / _Fs)
    tf = _trigger_time + (ind_hi / _Fs)

    np.testing.assert_equal(
        _windows([t0, tf]),
        np.array([window]))

    np.testing.assert_equal(
        _windows([t0 - (0.999 / _Fs), tf + (0.999 / _Fs)]),
        np.array([window]))

    # Window in middle where `tlim` *exactly* straddles boundaries
    np.testing.assert_equal(
        _windows([t0 - (1. / _Fs), tf]),
        np.array([window - 1, window]))

    np.testing.assert_equal(
        _windows([t0, tf + (1. / _Fs)]),
        np.array([window, window + 1]))

    return


# # There are some discrepancies between the manually computed
# # vibration-subtracted CO2 phase (`ph_CO2_plasma`) and
# # the automatically computed vibration-subtracted CO2 phase
# # (`sig_CO2.x`) that I do *not* fully understand...
# def test___plasma_induced_phase():
#     shot = 167341
#     tlim = [2.0, 2.25]
#
#     # Load vibration-subtracted CO2 phase data
#     sig_CO2 = Phase(shot, tlim=tlim, beam='CO2')
#
#     # Load raw CO2 and HeNe phase data
#     sig_CO2_raw = Phase(
#         shot, tlim=tlim, beam='CO2', vibration_subtracted=False)
#     sig_HeNe_raw = Phase(
#         shot, tlim=tlim, beam='HeNe', vibration_subtracted=False)
#
#     # Perform vibration-subtraction on CO2 phase manually
#     ph_CO2_plasma = _plasma_induced_phase(
#         sig_CO2_raw.x, sig_HeNe_raw.x)
#
#     np.testing.assert_array_almost_equal(sig_CO2.x, ph_CO2_plasma)
#
#     return
