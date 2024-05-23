from __future__ import annotations

from ._core import __doc__, __version__, add, subtract, plot_sine_wave, plot_cosine_wave, plot_square_wave, plot_sawtooth_wave, plot_dft_sinus, plot_dft_cosinus,signal_show

__all__ = ["__doc__", "__version__", "add", "subtract", "plot_sine_wave", "plot_cosine_wave", "plot_square_wave", "plot_sawtooth_wave", "plot_dft_sinus", "plot_dft_cosinus","signal_show"]

from scipy.io import wavfile
import numpy as np

def wizualization(plik_we):

    rate, file = wavfile.read(plik_we)
    if len(file.shape) > 1:
        file = file[:, 0]
    time = np.arange(1000) / rate
    signal = file[:1000]
    signal_show(time.tolist(), signal.tolist(), "First 1000 samples")
