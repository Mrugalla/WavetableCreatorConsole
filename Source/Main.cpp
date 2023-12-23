#include <JuceHeader.h>
#include <functional>

#define oopsie(x) jassert(!(x))

namespace wtc
{
    static constexpr float pi = 3.14159265359f;
    static constexpr float tau = pi * 2.f;
    
    inline float pitchToFreq(float pitch, float basePitch = 69.f, float xen = 12.f, float masterTune = 440.f) noexcept
    {
		return std::pow(2.f, (pitch - basePitch) / xen) * masterTune;
	}

    inline float freqToPitch(float freq, float basePitch = 69.f, float xen = 12.f, float masterTune = 440.f) noexcept
    {
        return std::log2(freq / masterTune) * xen + basePitch;
    }

    inline float dbToGain(float db) noexcept
    {
		return std::pow(10.f, db * .05f);
	}

    inline float gainToDb(float gain) noexcept
    {
        return std::log10(gain) * 20.f;
    }

    inline int twoPowXInverse(int num) noexcept
    {
        auto x = 0;
        while (num > 1)
        {
			num >>= 1;
			++x;
		}
		return x;
    }

    template<typename Float>
    inline Float modDesmos(Float a, Float b) noexcept
    {
        while (a < static_cast<Float>(0))
            a += b;
        while (a >= b)
            a -= b;
        return a;
    }

    struct BiquadFilter
    {
        enum Type { LP, HP, BP, AP, NOTCH, PEAK, LS, HS };
        static constexpr double DefaultCutoffFc = 1000. / 44100.;
        static constexpr double DefaultResonance = .1;

        BiquadFilter() :
            type(Type::BP),
            bCoef{ 0., 0., 0. },
            aCoef{ 0., 0. },
            w{ 0., 0. },
            omega(DefaultCutoffFc),
            cosOmega(0.), sinOmega(0.),
            Q(DefaultResonance), alpha(0.), A(1.),
            a{ 0., 0., 0. },
            b{ 0., 0., 0. }
        {
        }

        void prepare() noexcept
        {
            update();
        }

        void operator()(double* samples, int numSamples) noexcept
        {
            for (auto s = 0; s < numSamples; ++s)
                processSample(samples[s]);
        }

        double operator()(double sample) noexcept
        {
            return processSample(sample);
        }

        /* fc = hz / sampleRate */
        void setCutoffFc(double fc) noexcept
        {
            omega = static_cast<double>(tau) * fc;
        }

        void setResonance(double q) noexcept
        {
            Q = q;
        }

        void setType(Type _type) noexcept
        {
            type = _type;
        }

        void update() noexcept
        {
            cosOmega = cos(omega);
            sinOmega = sin(omega);

            switch (type)
            {
            case LP:
                alpha = sinOmega / (2. * Q);
                b[0] = (1. - cosOmega) / 2.;
                b[1] = 1. - cosOmega;
                b[2] = b[0];
                a[0] = 1. + alpha;
                a[1] = -2. * cosOmega;
                a[2] = 1. - alpha;
                break;
            case HP:
                alpha = sinOmega / (2. * Q);
                b[0] = (1. + cosOmega) / 2.;
                b[1] = -(1. + cosOmega);
                b[2] = b[0];
                a[0] = 1. + alpha;
                a[1] = -2. * cosOmega;
                a[2] = 1. - alpha;
                break;
            case BP:
                alpha = sinOmega * sinh(log(2.) / 2. * Q * omega / sinOmega);
                b[0] = sinOmega / 2.;
                b[1] = 0.;
                b[2] = -b[0];
                a[0] = 1. + alpha;
                a[1] = -2. * cosOmega;
                a[2] = 1. - alpha;
                break;
            case AP:
                alpha = sinOmega / (2. * Q);
                b[0] = 1. - alpha;
                b[1] = -2. * cosOmega;
                b[2] = 1. + alpha;
                a[0] = b[2];
                a[1] = b[1];
                a[2] = b[0];
                break;
            case NOTCH:
                alpha = sinOmega * sinh(log(2.) / 2. * Q * omega / sinOmega);
                b[0] = 1.;
                b[1] = -2. * cosOmega;
                b[2] = 1.;
                a[0] = 1. + alpha;
                a[1] = b[1];
                a[2] = 1. - alpha;
                break;
            case PEAK:
                alpha = sinOmega * sinh(log(2.) / 2. * Q * omega / sinOmega);
                b[0] = 1. + (alpha * A);
                b[1] = -2. * cosOmega;
                b[2] = 1. - (alpha * A);
                a[0] = 1. + (alpha / A);
                a[1] = b[1];
                a[2] = 1. - (alpha / A);
                break;
            case LS:
                alpha = sinOmega / 2. * sqrt((A + 1. / A) * (1. / Q - 1.) + 2.);
                b[0] = A * ((A + 1.) - ((A - 1.) * cosOmega) + (2. * sqrt(A) * alpha));
                b[1] = 2. * A * ((A - 1.) - ((A + 1.) * cosOmega));
                b[2] = A * ((A + 1.) - ((A - 1.) * cosOmega) - (2. * sqrt(A) * alpha));
                a[0] = ((A + 1.) + ((A - 1.) * cosOmega) + (2. * sqrt(A) * alpha));
                a[1] = -2. * ((A - 1.) + ((A + 1.) * cosOmega));
                a[2] = ((A + 1.) + ((A - 1.) * cosOmega) - (2. * sqrt(A) * alpha));
                break;
            case HS:
                alpha = sinOmega / 2. * sqrt((A + 1. / A) * (1. / Q - 1.) + 2.);
                b[0] = A * ((A + 1.) + ((A - 1.) * cosOmega) + (2. * sqrt(A) * alpha));
                b[1] = -2. * A * ((A - 1.) + ((A + 1.) * cosOmega));
                b[2] = A * ((A + 1.) + ((A - 1.) * cosOmega) - (2. * sqrt(A) * alpha));
                a[0] = ((A + 1.) - ((A - 1.) * cosOmega) + (2. * sqrt(A) * alpha));
                a[1] = 2. * ((A - 1.) - ((A + 1.) * cosOmega));
                a[2] = ((A + 1.) - ((A - 1.) * cosOmega) - (2. * sqrt(A) * alpha));
                break;
            }

            // Normalize filter coefficients
            const auto factor = 1. / a[0];

            aCoef[0] = a[1] * factor;
            aCoef[1] = a[2] * factor;

            bCoef[0] = b[0] * factor;
            bCoef[1] = b[1] * factor;
            bCoef[2] = b[2] * factor;
        }

        void copyFrom(const BiquadFilter& other) noexcept
        {
            type = other.type;
            for (auto i = 0; i < bCoef.size(); ++i)
                bCoef[i] = other.bCoef[i];
            for (auto i = 0; i < aCoef.size(); ++i)
                aCoef[i] = other.aCoef[i];
            omega = other.omega;
            cosOmega = other.cosOmega;
            sinOmega = other.sinOmega;
            Q = other.Q;
            alpha = other.alpha;
            A = other.A;
            for (auto i = 0; i < a.size(); ++i)
                a[i] = other.a[i];
            for (auto i = 0; i < b.size(); ++i)
                b[i] = other.b[i];
        }

        Type type;

        std::array<double, 3> bCoef; // b0, b1, b2
        std::array<double, 2> aCoef; // a1, a2
        std::array<double, 2> w; // delays

        double omega, cosOmega, sinOmega, Q, alpha, A;
        std::array<double, 3> a;
        std::array<double, 3> b;

        // DF-II impl
        double processSample(double sample) noexcept
        {
            auto y = bCoef[0] * sample + w[0];
            w[0] = bCoef[1] * sample - aCoef[0] * y + w[1];
            w[1] = bCoef[2] * sample - aCoef[1] * y;
            return y;
        }
    };

    using Buffer = juce::AudioBuffer<float>;

    struct Specs
    {
        Specs(int _tableSize = 2048, int _numTables = 256) :
            tableSize(_tableSize),
            numTables(_numTables),
            fftOrder(twoPowXInverse(tableSize)),
            fftSize(1 << fftOrder),
            fftSize2(fftSize * 2),
            tableSizeF(static_cast<float>(tableSize)),
            numTablesF(static_cast<float>(numTables)),
            fftSizeF(static_cast<float>(fftSize)),
            fftSize2F(static_cast<float>(fftSize2)),
            bins(),
            fifo(),
            fft(fftOrder),
            randomizer()
        {
            bins.resize(fftSize2, 0.f);
            fifo.resize(fftSize2, 0.f);
        }

        void clearFifo() noexcept
        {
            for (auto& f : fifo)
                f = std::complex<float>(0.f, 0.f);
        }

        void setMagAndPhase(float mag, float phase, int idx) noexcept
        {
            fifo[idx] = std::complex<float>
                (
                    mag * std::cos(phase),
                    mag * std::sin(phase)
                );
        }

        void addMagAndPhase(float mag, float phase, int idx) noexcept
        {
            fifo[idx] = std::complex<float>
                (
                    fifo[idx].real() + mag * std::cos(phase),
                    fifo[idx].imag() + mag * std::sin(phase)
                );
        }

        void randomizePhases(float amount = 1.f) noexcept
        {
			for (auto& f : fifo)
                if(f.real() != 0.f)
				    f = std::complex<float>(f.real(), rand() * tau * amount);
		}

        void makeSomeNoise() noexcept
        {
			for (auto i = 0; i < fftSize; ++i)
                fifo[i] = std::complex<float>(1.f, rand() * tau);
		}

        void applyWindowHann() noexcept
        {
            for (auto i = 0; i < fftSize; ++i)
            {
                const auto hann = .5f * (1.f - std::cos(tau * static_cast<float>(i) / fftSizeF));
                fifo[i] *= hann;
            }
        }

        void addFunction(float x, const std::function<float(float)>& magFunc, 
            const std::function<float(float)>& phaseFunc = [](float) { return 0.f; },
            const std::function<float(float)>& gainFunc = [](float) { return 1.f; }) noexcept
        {
            const auto pitchStart = freqToPitch(1.f);
            const auto pitchEnd = freqToPitch(tableSizeF * .5f);
            const auto pitch = pitchStart + magFunc(x) * (pitchEnd - pitchStart);
            const auto freq = pitchToFreq(pitch);
            const auto fc = freq / tableSizeF;
            const auto idx = static_cast<int>(std::round(fc * fftSizeF));
            if(idx < 0 || idx >= fftSize)
				return;
            addMagAndPhase(gainFunc(x), phaseFunc(x), idx);
        }

        void addFuncHarmonic(float x = 0.f) noexcept
        {
			addFunction
            (
                0.f,
				[x](float) { return x; }
			);
		}

        void normalizeFifo() noexcept
        {
            auto max = fifo[0].real();
            for (auto i = 1; i < fftSize2; ++i)
            {
                const auto a = std::abs(fifo[i].real());
                if (max < a)
                    max = a;
            }
            if (max != 0.f)
            {
                auto gain = 1.f / max;
                for (auto& f : fifo)
                    f *= gain;
            }
        }

        void performInverseFFT() noexcept
        {
            normalizeFifo();
            fft.perform(fifo.data(), bins.data(), true);
        }

        float getSpectralSound(int sampleIdx) const noexcept
        {
            return bins[sampleIdx].real();
        }

        float rand()
        {
            return randomizer.nextFloat();
        }

        int tableSize, numTables, fftOrder, fftSize, fftSize2;
        float tableSizeF, numTablesF, fftSizeF, fftSize2F;
        std::vector<std::complex<float>> bins;
        std::vector<std::complex<float>> fifo;
    private:
        juce::dsp::FFT fft;
        juce::Random randomizer;
    };

    namespace fftFX
    {
        inline void distort(Specs& specs, float gain, float distance = 2.f, bool preventAliasing = true) noexcept
        {
            for (auto i = static_cast<int>(specs.fftSizeF / distance); i >= 0; --i)
            {
                const auto fc = static_cast<float>(i) / specs.fftSize2F;
                const auto fc2 = fc * distance;
                if (fc2 >= 0.f && fc2 < .5f)
                {
                    if (!preventAliasing || fc2 < specs.fftSizeF - 1.f)
                    {
                        const auto dFloor = fc2 * specs.fftSize2F;
                        const auto d0 = static_cast<int>(dFloor);
                        const auto d1 = d0 + 1;
                        const auto dX = dFloor - static_cast<float>(d0);
                        const auto dMag = specs.fifo[i].real() * gain;
                        const auto dMag0 = dMag * (1.f - dX);
                        const auto dMag1 = dMag * dX;
                        const auto dPhs = specs.fifo[i].imag() * gain;
                        const auto dPhs0 = dPhs * (1.f - dX);
                        const auto dPhs1 = dPhs * dX;
                        specs.addMagAndPhase(dMag0, dPhs0, d0);
                        specs.addMagAndPhase(dMag1, dPhs1, d1);
                    }
                }
            }
        }

        inline void blurBrk(Specs& specs, float fcWidth)
        {
            auto nFifo = specs.fifo;

            for (auto i = 0; i < specs.fftSize; ++i)
            {
                const auto mag = specs.fifo[i].real() * .5f;
                const auto phs = specs.fifo[i].imag() * .5f;
                const auto iF = static_cast<float>(i);
                const auto fc = iF / specs.fftSize2F;
                const auto fc0 = fc - fcWidth;
                if (fc0 > 0.f)
                {
                    const auto idx0 = static_cast<int>(std::round(fc0 * specs.fftSize2F));
                    nFifo[idx0] += { mag, phs };
                }
                const auto fc1 = fc + fcWidth;
                if (fc1 < 1.f)
                {
                    const auto idx1 = static_cast<int>(std::round(fc1 * specs.fftSize2F));
                    nFifo[idx1] += { mag, phs };
                }
            }

            specs.fifo = nFifo;
        }

        inline void blur(Specs& specs,
            float fc, float reso, int slope) noexcept
        {
            BiquadFilter lpMag, lpPhase;
            lpMag.setType(BiquadFilter::Type::LP);
            lpMag.setCutoffFc(fc);
            lpMag.setResonance(reso);
            lpMag.update();
            lpPhase.setType(BiquadFilter::Type::LP);
            lpPhase.setCutoffFc(fc);
            lpPhase.setResonance(reso);
            lpPhase.update();

            auto dest = specs.fifo;
            auto src = specs.fifo.data();

            const auto lp = [&](int i)
            {
                auto mag = static_cast<double>(src[i].real());
                auto phs = static_cast<double>(src[i].imag());
                mag = lpMag(mag);
                phs = lpPhase(phs);
                dest[i] = { static_cast<float>(mag), static_cast<float>(phs) };
			};

            for (auto i = specs.fftSize - 1; i != 0; --i)
                lp(i);
            for (auto i = 0; i < specs.fftSize; ++i)
                lp(i);

            for (auto i = 1; i < slope; ++i)
            {
                for (auto j = specs.fftSize - 1; j != 0; --j)
                    lp(i);
                for (auto j = 0; j < specs.fftSize; ++j)
                    lp(i);
            }
        }
    }

    using PrepareFunc = std::function<void(Specs&, int/*wavetable index*/)>;
    using TableFunc = std::function<float(Specs& specs, float/*x*/, int/*sample index*/, int/*wavetable index*/)>;

    struct Creator
    {
        Creator(Specs& specs, const TableFunc& func, const PrepareFunc& prepare = [](Specs&, int){}) :
            tableSizeInv(1.f / static_cast<float>(specs.tableSize)),
            tableSize(specs.tableSize),
            numTables(specs.numTables),
            buffer(1, tableSize * numTables)
        {
            buffer.clear();
            auto samples = buffer.getWritePointer(0);

            for (auto i = 0; i < numTables; ++i)
            {
                prepare(specs, i);

                const auto start = i * tableSize;
                const auto end = (i + 1) * tableSize;
                for (auto s = start; s < end; ++s)
                {
                    const auto s0 = s % tableSize;
                    const auto sF = static_cast<float>(s0);
                    const auto x = 2.f * sF * tableSizeInv - 1.f;

                    auto sample = func(specs, x, s0, i);

                    if (std::isnan(sample) || std::isinf(sample))
                        sample = 0.f;

                    samples[s] = sample;
                }
            }
        }

        // POST PROCESSING (SINGLE)
        void dcOffsetEachTable() noexcept
        {
            auto samples = buffer.getWritePointer(0);

            for (auto i = 0; i < numTables; ++i)
            {
                const auto start = tableSize * i;
                const auto end = tableSize * (i + 1);

                auto average = 0.f;

                for (auto s = start; s < end; ++s)
                    average += samples[s];

                const auto offset = -average / static_cast<float>(end - start);

                for (auto s = start; s < end; ++s)
                    samples[s] += offset;
            }
        }

        void normalizeEachTable() noexcept
        {
            auto samples = buffer.getWritePointer(0);

            for (auto i = 0; i < numTables; ++i)
            {
                const auto start = tableSize * i;
                const auto end = tableSize * (i + 1);

                auto max = 0.f;

                for (auto s = start; s < end; ++s)
                {
                    const auto a = std::abs(samples[s]);
                    if (max < a)
                        max = a;
                }

                const auto gain = 1.f / max;

                for (auto s = start; s < end; ++s)
                    samples[s] *= gain;
            }
        }

        void applyHannWindowEachTable() noexcept
        {
            auto samples = buffer.getWritePointer(0);

            for (auto i = 0; i < numTables; ++i)
            {
                auto k = i * tableSize;

                for (auto j = 0; j < tableSize; ++j)
                {
					const auto hann = .5f * (1.f - std::cos(tau * static_cast<float>(j) * tableSizeInv));
					samples[k + j] *= hann;
				}
            }
        }

        void xFadeTablesSingle(float lenRel)
        {
            std::vector<float> table2;
            auto fadeLen = static_cast<int>(lenRel * static_cast<float>(tableSize));
            table2.resize(fadeLen);

            for (auto i = 0; i < numTables; ++i)
            {
                auto table = &buffer.getWritePointer(0)[i * tableSize];

                for (auto j = 0; j < fadeLen; ++j)
                    table2[j] = table[fadeLen - j - 1];

                for (auto j = 0; j < fadeLen; ++j)
                {
                    auto idx = tableSize - fadeLen + j;
                    auto x = (float)j * tableSizeInv;
                    table[idx] += x * (table2[j] - table[idx]);
                }
            }
        }

        // POST PROCESSING (MULTI)

        void xFadeTablesMulti() noexcept
        {
            for (auto i = 0; i < numTables - 1; ++i)
            {
                auto table0 = &buffer.getWritePointer(0)[i * tableSize];
                auto table1 = &buffer.getWritePointer(0)[(i + 1) * tableSize];

                for (auto j = 0; j < tableSize; ++j)
                {
                    auto x = static_cast<float>(j) * tableSizeInv;
                    table0[j] = table0[j] + x * (table1[j] - table0[j]);
                }
            }
        }

        void upsampleNumTables(int factor)
        {
            auto newNumTables = numTables * factor;
            Buffer newBuffer(1, tableSize * newNumTables);
            auto newSamples = newBuffer.getWritePointer(0);

            const auto samplesOld = buffer.getReadPointer(0);
            auto samplesNew = newBuffer.getWritePointer(0);

            for (auto i = 0; i < newNumTables; ++i)
            {
                auto iNew = i * tableSize;
                auto iOld = iNew / factor;
                auto iOld0 = iOld / tableSize * tableSize;
                auto x = static_cast<float>(iOld) * tableSizeInv;
                while (x >= 1.f)
                    --x;
                for (auto j = 0; j < tableSize; ++j)
                {
                    auto jNew = iNew + j;
                    auto jOld0 = iOld0 + j;
                    auto jOld1 = jOld0 + tableSize;
                    if (jOld1 >= buffer.getNumSamples())
                        jOld1 -= tableSize;
                    samplesNew[jNew] = samplesOld[jOld0] + x * (samplesOld[jOld1] - samplesOld[jOld0]);
                }
            }

            numTables = newNumTables;
            buffer = newBuffer;
        }

        void upsampleNumTablesTo(int newNumTables)
        {
            if (newNumTables == numTables)
                return;
			const auto factor = static_cast<int>(static_cast<float>(newNumTables) / static_cast<float>(numTables));
			upsampleNumTables(factor);
		}

        // INTERPOLATION NumTables
        // ...

        // EXPORT
        void makeWav(const juce::String& name)
        {
            const auto numSamples = buffer.getNumSamples();
            auto samples = buffer.getWritePointer(0);
            for (auto s = 0; s < numSamples; ++s)
                samples[s] = juce::jlimit(-1.f, 1.f, samples[s]);

            const auto path = juce::File::getSpecialLocation
            (
                juce::File::SpecialLocationType::userDesktopDirectory
            ).getFullPathName() + "\\WaveTableCreator\\";
            juce::File file(path);
            file.createDirectory();
			std::cout << "Created Folder: " << path << std::endl;
            const auto fullName = name + ".wav";
            file = juce::File(path + "\\" + fullName);
            if (file.existsAsFile())
                file.deleteFile();
            std::cout << "Created File: " << file.getFullPathName() << std::endl;
            juce::WavAudioFormat format;
            std::unique_ptr<juce::AudioFormatWriter> writer;
            writer.reset(format.createWriterFor(new juce::FileOutputStream(file),
                48000,
                1,
                32,
                {},
                0
            ));
            if (writer != nullptr)
            {
                writer->writeFromAudioSampleBuffer(buffer, 0, numSamples);
                std::cout << "Success :)" << std::endl;

                juce::File txtFile(path + "\\" + "FolderInfo.txt");
                if (txtFile.existsAsFile())
					txtFile.deleteFile();
                txtFile.create();
                txtFile.appendText("[" + juce::String(tableSize) + "]");

                // open folder
                juce::URL url(path);
                url.launchInDefaultBrowser();
                url = file.getFullPathName();
                url.launchInDefaultBrowser();
            }
            else
            {
				std::cout << "Failed :(" << std::endl;
            }
        }

    protected:
        float tableSizeInv;
        int tableSize;
        int numTables;
        Buffer buffer;
    };
}

int main(int, char*)
{
    using namespace wtc;

    Specs specs(2048, 64);

    PrepareFunc prepare = [](Specs& specs, int tableIdx)
    {
        specs.clearFifo();
        const auto x = static_cast<float>(tableIdx) / specs.numTablesF;
        
        for (auto a = 0; a < 3; ++a)
        {
            specs.addFuncHarmonic(0.f);

            specs.addFunction(x, [a](float x)
                {
                    auto aa = (float)(1 << (a * 2));
                    auto xx = x * aa;
                    while (xx > 1.f)
                        xx -= 1.f;
                    return std::tan(xx * aa) / aa;
                });
            fftFX::blur(specs, .5f, 1.f, 3);
            fftFX::distort(specs, .75f, 1 << (a + 1), false);
        }
        
        specs.performInverseFFT();
    };

    TableFunc func = [](Specs& specs, float x, int sampleIdx, int tableIdx)
    {
        return specs.getSpectralSound(sampleIdx);
    };

    Creator creator(specs, func, prepare);

    creator.dcOffsetEachTable();
    //creator.xFadeTablesMulti();
    creator.normalizeEachTable();
    creator.upsampleNumTablesTo(256);
    creator.makeWav("xxx");

    return 0;
}

/*

wavetable ideas and stuff:

/////////////////////////////

// triangle
T(x,f) = 2asin(sin(fxpi))/pi

//either
T(x,f)sqrt(1-abs(x))
//or
sum(n=1,N=128,T(x,n)/n)

/////////////////////////////

LINEAR SWEEP FFT PREPARE:

specs.clearFifo();
auto norm = (float)tableIdx / (float)specs.numTables;
auto idx = static_cast<int>(norm * static_cast<float>(specs.fftSize));
specs.fifo[idx] = 1.f;
specs.performInverseFFT();

/////////////////////////////

LOG SWEEP FFT PREPARE:

specs.clearFifo();
auto norm = (float)tableIdx / (float)specs.numTables;
auto pitch = 48.f + norm * 64.f;
auto freq = pitchToFreq(pitch);
auto normFreq = freq / 44100.f;
auto idx = static_cast<int>(normFreq * static_cast<float>(specs.fftSize));
specs.fifo[idx] = 1.f;
specs.performInverseFFT();

/////////////////////////////

NOISY LOG SWEEP FFT PREPARE:

specs.clearFifo();
auto norm = (float)tableIdx / (float)specs.numTables;
auto pitch = 24.f + norm * 127.f;
auto freq = pitchToFreq(pitch);
auto normFreq = freq / 44100.f;
auto idx = static_cast<int>(normFreq * static_cast<float>(specs.fftSize));
specs.fifo[idx] = 1.f;
auto diffuse = 128;
for (auto i = 1; i < diffuse; ++i)
{
    auto iF = static_cast<float>(i);
    auto iX = iF / static_cast<float>(diffuse);
    iX = 1.f - std::pow(iX, 1.f / 64.f);
    auto i0 = idx + i;
    auto i1 = idx - i;
    if(i0 < specs.fftSize)
        specs.fifo[i0] = iX;
    if(i1 >= 0)
        specs.fifo[i1] = iX;
}

specs.performInverseFFT();

/////////////////////////////

specs.clearFifo();
auto norm = (float)tableIdx / (float)specs.numTables;
auto pitch = 24.f + norm * 127.f;
auto freq = pitchToFreq(pitch);
auto normFreq = freq / 44100.f;
auto idx = static_cast<int>(normFreq * static_cast<float>(specs.fftSize));
specs.fifo[idx] = 1.f;
auto diffuse = 8;
for (auto i = 1; i < diffuse; ++i)
{
    auto iF = static_cast<float>(i);
    auto iX = iF / static_cast<float>(diffuse);
    iX = 1.f - std::pow(iX, 1.f / 32.f);
    auto i0 = idx + i;
    auto i1 = idx - i;
    if (i0 < specs.fftSize)
        specs.fifo[i0] = iX;
    if (i1 >= 0)
        specs.fifo[i1] = iX;
}

specs.randomizePhases();
specs.performInverseFFT();

/////////////////////////////

URANUS RINGS FFT PREPARE:

specs.clearFifo();
auto x = static_cast<float>(tableIdx) / specs.numTablesF;
auto pitchSweepStart = 12.f;
auto pitchSweepLength = 64.f;
auto numHarmonics = 32;
auto pitchDistance = 5.f;
for (auto i = 0; i < numHarmonics; ++i)
{
    const auto iF = static_cast<float>(i);
    const auto iPitchDistance = pitchDistance * iF;
    auto pitch = pitchSweepStart + x * pitchSweepLength + iPitchDistance;
    while(pitch >= 136.f)
        pitch -= 136.f;
    const auto freq = pitchToFreq(pitch);
    const auto fc = freq / 44100.f;
    const auto idx = static_cast<int>(std::round(fc * specs.fftSizeF * 2.f));
    specs.fifo[idx] = 1.f;
}

specs.randomizePhases();
specs.performInverseFFT();

/////////////////////////////

LOWEST PITCH SINE WAVE FFT PREPARE:

specs.clearFifo();
const auto x = static_cast<float>(tableIdx) / specs.numTablesF;
const auto lowestPitch = -36.f;
const auto freq = pitchToFreq(lowestPitch);
const auto fc = freq / specs.tableSizeF;
const auto idx = static_cast<int>(std::round(fc * specs.fftSize2F));
specs.fifo[idx] = 1.f;
specs.performInverseFFT();

/////////////////////////////

SAW WAVE FFT PREPARE:

specs.clearFifo();
const auto x = static_cast<float>(tableIdx) / specs.numTablesF;
const auto numHarmonics = specs.fftSize;
for (auto i = 1; i <= numHarmonics; ++i)
{
    const auto iF = static_cast<float>(i);
    const auto freq = iF;
    const auto fc = freq / specs.tableSizeF;
    if (fc < .5f)
    {
        const auto idx = static_cast<int>(std::round(fc * specs.fftSize2F));
        auto mag = 1.f / iF;
        auto phase = 0.f;// i % 2 == 0 ? 0 : pi;
        specs.setMagAndPhase(mag, phase, idx);
    }
}
specs.performInverseFFT();

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

/////////////////////////////

TODO:

being able to create sub paths by typing it into the string of createWav
export creates helper txt file for serum

*/