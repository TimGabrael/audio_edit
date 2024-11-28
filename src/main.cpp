#include <iostream>
#include <complex>
#include "AudioFile.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <functional>


template<typename T> 
std::vector<std::complex<T>> CalculateDFT(T* data, size_t size) {
    std::vector<std::complex<T>> output(size);

    const T n = static_cast<T>(size);
    const T inv_n = 1.0 / n;

    const size_t half_size = size / 2 + 1;
    for(size_t i = 0; i < half_size; ++i) {
        T cur_real = 0.0;
        T cur_imag = 0.0;
        for(size_t j = 0; j < size; ++j) {
            T freq = 2.0 * M_PI * static_cast<T>(i) * static_cast<T>(j) * inv_n;
            cur_real += data[j] * std::cos(freq);
            cur_imag -= data[j] * std::sin(freq);
        }
        output.at(i).real(cur_real);
        output.at(i).imag(cur_imag);

        // use the symmetry of the 1d real dft
        if(i > 0 && i < (half_size - 1)) {
            output.at(output.size() - i) = std::conj(output.at(i));
        }
    }
    return output;
}

template<typename T> 
std::vector<std::complex<T>> CalculateInverseDFT(std::complex<T>* data, size_t size) {
    std::vector<std::complex<T>> output(size);
    const T n = static_cast<T>(size);
    const T inv_n = 1.0 / n;

    for(size_t i = 0; i < size; ++i) {
        T cur_real = 0.0;
        T cur_imag = 0.0;
        for(size_t j = 0; j < size; ++j) {
            T freq = 2.0 * M_PI * static_cast<T>(i) * static_cast<T>(j) * inv_n;
            T cos_val = std::cos(freq);
            T sin_val = std::sin(freq);
            cur_real += data[j].real() * cos_val - data[j].imag() * sin_val;
            cur_imag += data[j].real() * sin_val + data[j].imag() * cos_val;
        }
        output.at(i).real(cur_real * inv_n);
        output.at(i).imag(cur_imag * inv_n);
    }
    return output;
}

template<typename T>
T WrapPhase(T phase_in) {
    if(phase_in >= 0.0) {
        return std::fmod(phase_in + M_PI, 2.0 * M_PI) - M_PI;
    }
    return std::fmod(phase_in - M_PI, -2.0 * M_PI) + M_PI;
}

template<typename T, size_t fft_size, size_t hop_size>
struct PitchShiftData {
    T last_input_phases[fft_size] = {};
    T last_output_phases[fft_size] = {};

    T analysis_magnitudes[fft_size / 2 + 1] = {};
    T analysis_frequencies[fft_size / 2 + 1] = {};
    T synthesis_magnitudes[fft_size / 2 + 1] = {};
    T synthesis_frequencies[fft_size / 2 + 1] = {};

    double pitch_shift = 1.0;

    PitchShiftData(double pitch_shift) : pitch_shift(pitch_shift) {
    }
};

template<typename T, size_t fft_size, size_t hop_size>
void PitchShiftWorkerFunction(std::vector<std::complex<T>>& fft_data, PitchShiftData<T, fft_size, hop_size>& pitch_shift_data) {
    const T inv_fft_size = 1.0 / static_cast<T>(fft_data.size());
    const T inv_hop_2_pi = 1.0 / (hop_size * 2.0 * M_PI); 

    const size_t half_fft_size = fft_data.size() / 2 + 1;

    for(size_t i = 0; i < half_fft_size; ++i) {
        const T amplitude = std::abs(fft_data.at(i));
        const T phase = std::atan2(fft_data.at(i).imag(), fft_data.at(i).real());

        T phase_diff = phase - pitch_shift_data.last_input_phases[i];

        const T bin_centre_frequency = 2.0 * M_PI * static_cast<T>(i) * inv_fft_size;
        phase_diff = WrapPhase(phase_diff - bin_centre_frequency * hop_size);
        
        const T bin_deviation = phase_diff * static_cast<T>(fft_data.size()) * inv_hop_2_pi;

        pitch_shift_data.analysis_frequencies[i] = static_cast<T>(i) + bin_deviation;
        pitch_shift_data.analysis_magnitudes[i] = amplitude;

        pitch_shift_data.last_input_phases[i] = phase;
    }

    for(size_t i = 0; i < half_fft_size; ++i) {
        pitch_shift_data.synthesis_magnitudes[i] = 0.0;
        pitch_shift_data.synthesis_frequencies[i] = 0.0;
    }
    for(size_t i = 0; i < half_fft_size; ++i) {
        const size_t new_bin = static_cast<size_t>(std::floor(i * pitch_shift_data.pitch_shift + 0.5));
        if(new_bin < half_fft_size) {
            pitch_shift_data.synthesis_magnitudes[new_bin] += pitch_shift_data.analysis_magnitudes[i];
            pitch_shift_data.synthesis_frequencies[new_bin] = pitch_shift_data.analysis_frequencies[i] * pitch_shift_data.pitch_shift;
        }
    }

    for(size_t i = 0; i < half_fft_size; ++i) {
        const T amplitude = pitch_shift_data.synthesis_magnitudes[i];
        const T bin_deviation = pitch_shift_data.synthesis_frequencies[i] - static_cast<T>(i);

        T phase_diff = bin_deviation * 2.0 * M_PI * hop_size * inv_fft_size;

        const T bin_centre_frequency = 2.0 * M_PI * static_cast<T>(i) * inv_fft_size;

        phase_diff += bin_centre_frequency * hop_size;

        const T out_phase = WrapPhase(pitch_shift_data.last_output_phases[i] + phase_diff);
        fft_data.at(i) = std::polar(amplitude, out_phase);

        if(i > 0 && i < (half_fft_size - 1)) {
            fft_data.at(fft_data.size() - i) = std::conj(fft_data.at(i));
        }
        pitch_shift_data.last_output_phases[i] = out_phase;
    }


}

template<typename T>
void SetHannWindow(std::vector<T>& frame) {
    const T inv_n_1 = 1.0 / static_cast<T>(frame.size() - 1);
    for (size_t i = 0; i < frame.size(); ++i) {
        frame[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i * inv_n_1));
    }
}
template<typename T>
std::vector<T> PhaseVocoder(const std::vector<T>& input, std::function<void(std::vector<std::complex<T>>&)> func, size_t fft_size, size_t hop_size) {
    if(input.size() < fft_size) {
        std::cout << "[Error in PhaseVocoder]: Input is to small, has to be at least 1x fft_size" << std::endl;
        return {};
    }
    const size_t frame_count = (input.size() - fft_size) / hop_size + 1;

    std::vector<std::complex<T>> prev_phase(fft_size, {0.0f, 0.0f});
    std::vector<T> synth_buf(input.size(), 0.0);

    std::vector<T> window(fft_size, 1.0);
    SetHannWindow(window);

    for(size_t frame = 0; frame < frame_count; ++frame) {
        std::vector<T> input_frame;
        input_frame.reserve(fft_size);
        for(size_t i = 0; i < fft_size; ++i) {
            const size_t idx = i + frame * hop_size;
            if(idx < input.size()) {
                input_frame.push_back(input.at(idx) * window.at(i));
            }
        }
        if(input_frame.empty()) {
            break;
        }

        std::vector<std::complex<T>> freq_domain = CalculateDFT(input_frame.data(), input_frame.size());

        func(freq_domain);

        std::vector<std::complex<T>> time_domain = CalculateInverseDFT(freq_domain.data(), freq_domain.size());
        for(size_t i = 0; i < time_domain.size(); ++i) {
            const size_t idx = frame * hop_size + i;
            if(idx < synth_buf.size()) {
                synth_buf.at(idx) += time_domain.at(i).real() * window.at(i);
            }
        }
    }
    return synth_buf;
}

template<typename T, size_t fft_size, size_t hop_size>
std::vector<T> PitchShift(const std::vector<T>& input, double pitch_shift) {
    PitchShiftData<T, fft_size, hop_size> shift_data(pitch_shift);
    std::function<void(std::vector<std::complex<T>>&)> pitch_shift_worker = [&shift_data](std::vector<std::complex<T>>& v) {
        PitchShiftWorkerFunction<float, fft_size, hop_size>(v, shift_data);
    };
    std::vector<T> vocoded_data = PhaseVocoder<float>(input, pitch_shift_worker, fft_size, hop_size);
    return vocoded_data;
}

template<typename T>
std::vector<T> DefaultPitchShift(const std::vector<T>& input, double pitch_shift) {
    return PitchShift<T, 0x400, 0x100>(input, pitch_shift);
}

int main() {
    double pitch_shift = std::pow(2.0, 4.0 / 12.0); // 4 semi-tones up
    AudioFile<float> audio_data("assets/input.wav");
    audio_data.printSummary();

    if(audio_data.samples.size() > 0)  {
        std::vector<float> vocoded_data = DefaultPitchShift(audio_data.samples.at(0), pitch_shift);
        for(auto& samples : audio_data.samples) {
            samples.resize(vocoded_data.size());
            for(size_t i = 0; i < vocoded_data.size(); ++i) {
                samples.at(i) = vocoded_data.at(i);
            }
        }
        audio_data.save("assets/pitch_shifted.wav");
    }
    else {
        std::cout << "input has no samples" << std::endl;
        system("pause");
    }


    return 0;
}

