#include <JuceHeader.h>
#include <functional>

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

        void addFunction(const std::function<float(float)>& magFunc, 
            const std::function<float(float)>& phaseFunc, 
            const std::function<float(float)>& gainFunc, float x) noexcept
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

        void performInverseFFT() noexcept
        {
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

        // POST PROCESSING (MULTI)

        void xFadeTables() noexcept
        {
            for (auto i = 0; i < numTables - 1; ++i)
            {
                auto table0 = buffer.getWritePointer(0) + i * tableSize;
                auto table1 = buffer.getWritePointer(0) + (i + 1) * tableSize;

                for (auto j = 0; j < tableSize; ++j)
                {
                    auto x = static_cast<float>(j) / tableSize;
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

    Specs specs(2048, 256);

    PrepareFunc prepare = [](Specs& specs, int tableIdx)
    {
        specs.clearFifo();
        const auto x = static_cast<float>(tableIdx) / specs.numTablesF;
        specs.addFunction
        (
            [](float) { return 0.f; },
            [](float) { return 0.f; },
            [](float) { return 1.f; },
            x
        );
        auto len = 256.f;
        for(auto i = 1.f; i <= len; ++i)
            specs.addFunction
            (
                [i, len](float x)
                {
                    const auto X = i / len;
                    return -sin((1 + X) * cos((1 + 3 * X) * sin(1 + .0125 * X) * x * pi)) / 2.f + .5f;
                },
                [i, len](float x)
                {
                    const auto X = i / len;
                    return -cos((1 + X) * sin((1 + 3 * X) * cos(1 + .0125 * X) * x * pi)) / 2.f + .5f;
                },
                [i, len](float) { return 1.f; },
                x
            );
        specs.performInverseFFT();
    };

    TableFunc func = [](Specs& specs, float x, int sampleIdx, int tableIdx)
    {
        return specs.getSpectralSound(sampleIdx);
    };

    Creator creator(specs, func, prepare);

    //creator.upsampleNumTables(1);
    creator.dcOffsetEachTable();
    creator.normalizeEachTable();
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