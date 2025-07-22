import numpy as np
import matplotlib.pyplot as plt
from typing import Callable, List, Tuple, Literal, Union, Sequence

########################################################################################################################
BinSpec = Union[
    Tuple[float, float],            # (start, duration)
    Sequence[Tuple[float, float]]   # list/tuple of (start, duration) pairs
]

def waveform_generator(
    frequency: float,
    length: float,
    dt: float,
    envelope: Union[str, Callable[[np.ndarray], np.ndarray]] = "gaussian",
    phase: float = 0.0,
    Amplitude: float = 1.0,
    rise_time: float = 0.0,
    complex_output: bool = False,
    bins: BinSpec | None = None,
):
    """
    Generate an arbitrary waveform modulated by a specified envelope.
    
    Parameters
    ----------
    frequency : float (not a Angluar frequency)
        Carrier frequency (Hz).
    length : float
        Total duration of waveform (s).
    dt : float
        Time resolution (s).
    envelope : str | callable
        "gaussian", "square", "flat_top", or custom f(t)->env.
    phase : float, optional
        Carrier phase offset (rad).
    complex_output : bool, optional
        True -> exp(iωt); False -> cos(ωt).
    bins : (start, duration) or list thereof, optional
        Time windows to zero-out the envelope.  Example:
        bins = [(200e-9, 50e-9), (600e-9, 100e-9)]
    
    Returns
    -------
    t, waveform : np.ndarray, np.ndarray
    """
    t = np.arange(0, length+dt, dt)
    N = t.size

    # ─────────── envelope ───────────
    if isinstance(envelope, str):
        if envelope == "gaussian":
            """
            Gaussian envelope with standard deviation set to 1/6 of the total length.
            This ensures the envelope smoothly transitions to zero at the edges.
            """
            sigma  = length / 6
            center = length / 2
            env = np.exp(-0.5 * ((t - center) / sigma) ** 2)
        elif envelope == "square":
            """
            Square envelope that is 1 for the entire duration.
            """
            env = np.ones_like(t)
        elif envelope == "flat_top":
            """
            Flat-top envelope with a cosine rise/fall at the edges.
            The rise_time determines the width of the cosine transition,
            defaults to 1/10 of the total length if not specified.
            """
            if rise_time == 0:
                rise_time = length / 10
            rise_n      = max(int(np.round(rise_time / dt)), 1)
            edge        = 0.5 * (1 - np.cos(np.linspace(0, np.pi, rise_n)))
            env         = np.ones_like(t)
            env[:rise_n]  = edge
            env[-rise_n:] = edge[::-1]
        else:
            raise ValueError(f"Unsupported envelope: {envelope!r}")
    elif callable(envelope):
        """
        Custom envelope function that takes time array t and returns an envelope array.
        The function must return an array of the same length as t.
        """
        env = envelope(t)
        if env.size != N:
            raise ValueError("Custom envelope must match t in length.")
    else:
        raise TypeError("Envelope must be str or callable.")

    # ─────────── zero-amplitude bins ───────────
    if bins is not None:
        # accept a single tuple or an iterable of tuples
        bins_iter = bins if isinstance(bins, Sequence) and bins and isinstance(bins[0], tuple) else [bins]
        for start, dur in bins_iter:
            if start < 0 or dur < 0 or start + dur > length:
                raise ValueError(f"Invalid bin window (start={start}, dur={dur}).")
            env[(t >= start) & (t < start + dur)] = 0.0

    # ─────────── carrier ───────────
    phase_term = 2 * np.pi * frequency * t + phase
    carrier    = np.exp(1j * phase_term) if complex_output else np.sin(phase_term)

    # ─────────── Amplitude ───────────
    if Amplitude != 1.0:
        env = Amplitude * env

    waveform = env * carrier

    # ─────────── function form ───────────
    """
    this is for qutip.mesolve() compatibility.
    """

    def make_awg_hold(t_samples: np.ndarray, w_samples: np.ndarray):
        """
        Returns a function w_awg(t, args) that emulates
        an AWG zero–order hold on (t_samples, w_samples).
        """
        def w_awg(t_in, args=None):
            # find the index i such that t_samples[i] <= t_in < t_samples[i+1]
            idx = np.searchsorted(t_samples, t_in, side='right') - 1
            # clamp to valid range
            idx = max(0, min(idx, len(w_samples)-1))
            return w_samples[idx]
        return w_awg
    
    waveform_awg = make_awg_hold(t, waveform)

    return t, waveform, waveform_awg


    t_interp = np.arange(0, length+dt/2, dt/2)
    w_interp = np.interp(t_interp, t, waveform)


########################################################################################################################
CombineMode = Literal["concatenate", "sum"]
ReturnType   = Literal["numpy", "list"]

def combine_waveforms(
    waveform_list: List[Tuple[np.ndarray, np.ndarray]],
    mode: CombineMode = "concatenate",
    return_type: ReturnType = "numpy",
) -> Tuple[Union[np.ndarray, list], Union[np.ndarray, list]]:
    """
    Combine multiple (t, wave) tuples produced by `waveform_generator`.

    Parameters
    ----------
    waveform_list : list of (t, wave)
        Each entry is exactly what `waveform_generator` returns (NumPy arrays).
    mode : {'concatenate', 'sum'}, optional
        * 'concatenate' — stitches pulses end-to-end in time order.
        * 'sum'         — requires equal-length, equal-dt arrays and returns their element-wise sum.
    return_type : {'numpy', 'list'}, optional
        Choose 'list' to get Python lists back.

    Returns
    -------
    t_combined, wave_combined
        The unified time axis and waveform.
    """
    if not waveform_list:
        raise ValueError("waveform_list must not be empty.")

    if mode == "concatenate":
        # Stitch pulses sequentially
        t_segments, w_segments = zip(*waveform_list)
        # Ensure constant dt within each segment
        dt = np.diff(t_segments[0]).mean()
        # Re-index each segment so time is continuous
        t_offsets = [0.0]
        for seg in t_segments[:-1]:
            t_offsets.append(t_offsets[-1] + seg[-1] + dt)
        t_stitched = [t + off for t, off in zip(t_segments, t_offsets)]
        t_combined = np.concatenate(t_stitched)
        w_combined = np.concatenate(w_segments)

    elif mode == "sum":
        # Overlay pulses; all lengths & grids must match
        t_refs, w_refs = zip(*waveform_list)
        for t in t_refs[1:]:
            if t.size != t_refs[0].size or not np.allclose(t, t_refs[0]):
                raise ValueError("All waveforms must share identical time axes for 'sum' mode.")
        t_combined = t_refs[0]
        w_combined = np.sum(w_refs, axis=0)

    else:
        raise ValueError("mode must be 'concatenate' or 'sum'.")

    if return_type == "list":
        return t_combined.tolist(), w_combined.tolist()
    elif return_type == "numpy":
        return t_combined, w_combined
    else:
        raise ValueError("return_type must be 'numpy' or 'list'.")