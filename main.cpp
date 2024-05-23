#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <pybind11/stl.h>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <iomanip>
#include <pybind11/numpy.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

void signal_show(const std::vector<double>& time, const std::vector<double>& signal, const std::string& name = "") {
    matplot::plot(time, signal);
    if (!name.empty()) {
        matplot::title(name);
    }
    matplot::xlabel("time [s]");
    matplot::show();
}

void add_noise(std::vector<double>& values) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 0.1); // staly poziom szumu 0.1
    for (size_t i = 0; i < values.size(); ++i) {
        values[i] += distribution(gen);
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
    for (int i = 0; i < dft_result.size(); ++i) {
        magnitudes.push_back(std::abs(dft_result[i]));
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
       
        Functions:
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
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
    )pbdoc");

    m.def("plot_sine_wave", &plot_sine_wave,R"pbdoc(
        Plotting sine wave

        Parameters:
        frequency : double
        amplitude : double
        num_points : int
        noisy : bool - adding noise to signal

    )pbdoc",py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_cosine_wave", &plot_cosine_wave, R"pbdoc(
        Plotting cosine wave

        Parameters:
        frequency : double
        amplitude : double
        num_points : int
        noisy : bool - adding noise to signal

    )pbdoc", py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_square_wave", &plot_square_wave, R"pbdoc(
        Plotting square wave

        Parameters:
        frequency : double
        amplitude : double
        num_points : int
        noisy : bool - adding noise to signal

    )pbdoc", py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_sawtooth_wave", &plot_sawtooth_wave, R"pbdoc(
        Plotting sawtooth wave

        Parameters:
        frequency : double
        amplitude : double
        num_points : int
        noisy : bool - adding noise to signal

    )pbdoc", py::arg("frequency"), py::arg("amplitude"), py::arg("num_points"), py::arg("noisy"));

    m.def("plot_dft_sinus", &plot_dft_sinus, R"pbdoc(
        Plot the discrete fourier transform (DFT) of sinus

        Parameters:
        frequency : double
        num_points : int

    )pbdoc", py::arg("frequency"), py::arg("num_points"));

    m.def("plot_dft_cosinus", &plot_dft_cosinus, R"pbdoc(
        Plot the discrete fourier transform (DFT) of a cosinus

        Parameters:
        frequency : double
        num_points : int

    )pbdoc", py::arg("frequency"), py::arg("num_points"));

   m.def("signal_show", &signal_show, py::arg("time"), py::arg("signal"), py::arg("title") = "");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}

