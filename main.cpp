#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <pybind11/stl.h>
#include <cmath>
#include <vector>
#include <complex>
#include <random>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

void add_noise(std::vector<double>& values) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 0.1); // staly poziom szumu 0.1
    for (auto& value : values) {
        value += distribution(gen);
    }
}

std::vector<std::complex<double>> calculate_dft(const std::vector<double>& signal) {
    int N = signal.size();
    std::vector<std::complex<double>> dft(N, {0, 0});

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum = {0, 0};
        for (int n = 0; n < N; ++n) {
            double angle = 2 * M_PI * k * n / N;
            sum += signal[n] * std::exp(std::complex<double>(0, -angle));
        }
        dft[k] = sum;
    }

    return dft;
}

void plot_dft(const std::vector<double>& values) {
    std::vector<std::complex<double>> dft_result = calculate_dft(values);
    std::vector<double> magnitudes;
    for (const auto& value : dft_result) {
        magnitudes.push_back(std::abs(value));
    }
    matplot::plot(magnitudes);
    matplot::show();
}

void plot_dft_sinus(double frequency, int num_points)
{
    double amplitude = 1.0;
    std::vector<double> values(num_points);
    for (int i = 0; i < num_points; ++i) {
        values[i] = amplitude * std::sin(2 * M_PI * frequency * i / num_points);
    }
    
    plot_dft(values);
}

void plot_dft_cosinus(double frequency, int num_points)
{
    double amplitude = 1.0;
    std::vector<double> values(num_points);
    for (int i = 0; i < num_points; ++i) {
        values[i] = amplitude * std::cos(2 * M_PI * frequency * i / num_points);
    }
    
    plot_dft(values);
}

void plot_sine_wave(double frequency, double amplitude, int num_points, bool noisy)
{
    std::vector<double> values(num_points);
    for (int i = 0; i < num_points; ++i) {
        values[i] = amplitude * std::sin(2 * M_PI * frequency * i / num_points);
    }
    if (noisy) {
        add_noise(values);
    }
    matplot::plot(values);
    matplot::show();
}

void plot_cosine_wave(double frequency, double amplitude, int num_points, bool noisy)
{
    std::vector<double> values(num_points);
    for (int i = 0; i < num_points; ++i) {
        values[i] = amplitude * std::cos(2 * M_PI * frequency * i / num_points);
    }
    if (noisy) {
        add_noise(values);
    }
    matplot::plot(values);
    matplot::show();
}

void plot_square_wave(double frequency, double amplitude, int num_points, bool noisy)
{
    std::vector<double> values(num_points);
    for (int i = 0; i < num_points; ++i) {
        values[i] = (std::sin(2 * M_PI * frequency * i / num_points) >= 0 ? amplitude : -amplitude);
    }
    if (noisy) {
        add_noise(values);
    }
    matplot::plot(values);
    matplot::show();
}

void plot_sawtooth_wave(double frequency, double amplitude, int num_points, bool noisy)
{
    std::vector<double> values(num_points);
    for (int i = 0; i < num_points; ++i) {
        values[i] = (2 * (i * frequency / num_points - std::floor(0.5 + i * frequency / num_points)) * amplitude);
    }
    if (noisy) {
        add_noise(values);
    }
    matplot::plot(values);
    matplot::show();
}

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
       Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
           plot_sine_wave
           plot_cosine_wave
           plot_square_wave
           plot_sawtooth_wave
           plot_dft_sinus
           plot_dft_cosinus

    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("plot_sine_wave", &plot_sine_wave, R"pbdoc(
        Generate and plot a sine wave

        Parameters
        ----------
        frequency : double
            Frequency of the sine wave
        amplitude : double
            Amplitude of the sine wave
        num_points : int
            Number of points in the wave
        noisy : bool
            Whether to add noise to the wave
    )pbdoc", py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_cosine_wave", &plot_cosine_wave, R"pbdoc(
        Generate and plot a cosine wave

        Parameters
        ----------
        frequency : double
            Frequency of the cosine wave
        amplitude : double
            Amplitude of the cosine wave
        num_points : int
            Number of points in the wave
        noisy : bool
            Whether to add noise to the wave
    )pbdoc", py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_square_wave", &plot_square_wave, R"pbdoc(
        Generate and plot a square wave

        Parameters
        ----------
        frequency : double
            Frequency of the square wave
        amplitude : double
            Amplitude of the square wave
        num_points : int
            Number of points in the wave
        noisy : bool
            Whether to add noise to the wave
    )pbdoc", py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_sawtooth_wave", &plot_sawtooth_wave, R"pbdoc(
        Generate and plot a sawtooth wave

        Parameters
        ----------
        frequency : double
            Frequency of the sawtooth wave
        amplitude : double
            Amplitude of the sawtooth wave
        num_points : int
            Number of points in the wave
        noisy : bool
            Whether to add noise to the wave
    )pbdoc", py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_dft_sinus", &plot_dft_sinus, R"pbdoc(
        Plot the discrete Fourier transform (DFT) of a sinusoidal wave

        Parameters
        ----------
        frequency : double
            Frequency of the sinusoidal wave
        num_points : int
            Number of points in the wave
    )pbdoc", py::arg("frequency"), py::arg("num_points"));

    m.def("plot_dft_cosinus", &plot_dft_cosinus, R"pbdoc(
        Plot the discrete Fourier transform (DFT) of a cosinusoidal wave

        Parameters
        ----------
        frequency : double
            Frequency of the cosinusoidal wave
        num_points : int
            Number of points in the wave
    )pbdoc", py::arg("frequency"), py::arg("num_points"));

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

