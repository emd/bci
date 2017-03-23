from nose import tools
import numpy as np
from bci.signal import (
    _closest_digitized_point, _crop, _windows,
    trigger_time, Fs, Npts_per_window, Nwindows, Npts_total)


def test__closest_digitized_point():
    # Times less than `trigger_time` should point to first digitized point
    tools.assert_equal(
        _closest_digitized_point(trigger_time - 1),
        0)

    # Time at `trigger_time` should point to first digitized point
    tools.assert_equal(
        _closest_digitized_point(trigger_time),
        0)

    # Time at last time should point to last digitized point
    tf = trigger_time + (Npts_total / Fs)
    tools.assert_equal(
        _closest_digitized_point(tf),
         Npts_total - 1)

    # Times after last time should point to last digitized point
    tools.assert_equal(
        _closest_digitized_point(tf + 1),
         Npts_total - 1)

    # Somewhere in the middle
    divisor = 2
    t = trigger_time + ((Npts_total // divisor) * (1. / Fs))
    tools.assert_equal(
        _closest_digitized_point(t),
        (Npts_total // divisor))

    # Check that things work when time does *not* correspond
    # to a digitizer time stamp
    eps = np.finfo('float64').eps

    tools.assert_equal(
        _closest_digitized_point(t + (0.5 / Fs)),  # no change
        (Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t + ((0.5 + eps) / Fs)),  # rounds up
        (Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t - ((0.5 - eps) / Fs)),  # no change
        (Npts_total // divisor))

    tools.assert_equal(
        _closest_digitized_point(t - (0.5 / Fs)),  # rounds down
        (Npts_total // divisor))

    return


# def test__crop():
#     pass
# 
# 
def test__windows():
    # Lowest window
    tlim = [trigger_time, trigger_time + (1000. / Fs)]
    np.testing.assert_equal(
        _windows(tlim),
        np.array([0]))

    # Highest window
    tlim = np.array([Npts_total - 1000, Npts_total]) * (1. / Fs)
    np.testing.assert_equal(
        _windows(tlim),
        np.array([Nwindows - 1]))

    # Full set of windows
    tlim = trigger_time + (np.array([0, Npts_total]) * (1. / Fs))
    np.testing.assert_equal(
        _windows(tlim),
        np.arange(Nwindows))

    # Window in middle where `tlim` *exactly* sits on boundaries
    window = 3
    ind_lo = window * Npts_per_window              # Lower bound, global ind
    ind_hi = ((window + 1) * Npts_per_window) - 1  # Upper bound, global ind

    t0 = trigger_time + (ind_lo / Fs)
    tf = trigger_time + (ind_hi / Fs)

    np.testing.assert_equal(
        _windows([t0, tf]),
        np.array([window]))

    # Window in middle where `tlim` *exactly* straddles boundaries
    np.testing.assert_equal(
        _windows([t0 - (1. / Fs), tf]),
        np.array([window - 1, window]))

    np.testing.assert_equal(
        _windows([t0, tf + (1. / Fs)]),
        np.array([window, window + 1]))

    return
