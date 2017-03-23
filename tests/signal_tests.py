from nose import tools
import numpy as np
from bci.signal import (
    _closest_digitized_point, _crop, _windows,
    _trigger_time, _Fs, _Npts_per_window, _Nwindows, _Npts_total)


def test__closest_digitized_point():
    # Times less than `_trigger_time` should point to first digitized point
    tools.assert_equal(
        _closest_digitized_point(_trigger_time - 1),
        0)

    # Time at `_trigger_time` should point to first digitized point
    tools.assert_equal(
        _closest_digitized_point(_trigger_time),
        0)

    # Time at last time should point to last digitized point
    tf = _trigger_time + (_Npts_total / _Fs)
    tools.assert_equal(
        _closest_digitized_point(tf),
         _Npts_total - 1)

    # Times after last time should point to last digitized point
    tools.assert_equal(
        _closest_digitized_point(tf + 1),
         _Npts_total - 1)

    # Somewhere in the middle
    divisor = 2
    t = _trigger_time + ((_Npts_total // divisor) * (1. / _Fs))
    tools.assert_equal(
        _closest_digitized_point(t),
        (_Npts_total // divisor))

    # Check that things work when time does *not* correspond
    # to a digitizer time stamp
    eps = np.finfo('float64').eps

    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / _Fs)),  # no change
        (_Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t + ((0.5 + eps) / _Fs)),  # rounds up
        (_Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t - ((0.5 - eps) / _Fs)),  # no change
        (_Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t - (0.5 / _Fs)),  # rounds down
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
    t0 = _trigger_time + (gstart / _Fs)
    tf = _trigger_time + (gstop / _Fs)
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

    # Window in middle where `tlim` *exactly* straddles boundaries
    np.testing.assert_equal(
        _windows([t0 - (1. / _Fs), tf]),
        np.array([window - 1, window]))

    np.testing.assert_equal(
        _windows([t0, tf + (1. / _Fs)]),
        np.array([window, window + 1]))

    return
