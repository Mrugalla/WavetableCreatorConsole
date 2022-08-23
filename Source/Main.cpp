#include <JuceHeader.h>
#include <functional>

namespace wtc
{
    static constexpr float pi = 3.14159265359f;
    static constexpr float tau = pi * 2.f;
    
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

    using TableFunc = std::function<float(float/*x*/, int/*wavetable index*/)>;

    struct Creator
    {
        Creator(int _tableSize, int _numTables, const TableFunc& func) :
            tableSizeInv(1.f / static_cast<float>(_tableSize)),
            tableSize(_tableSize),
            numTables(_numTables),
            buffer(1, tableSize * numTables)
        {
            buffer.clear();

            for (auto i = 0; i < numTables; ++i)
            {
                const auto start = i * tableSize;
                const auto end = (i + 1) * tableSize;

                auto samples = buffer.getWritePointer(0);

                for (auto s = start; s < end; ++s)
                {
                    const auto s0 = s % tableSize;
                    const auto sF = static_cast<float>(s0);
                    const auto x = 2.f * sF * tableSizeInv - 1.f;

                    auto sample = func(x, i);

                    if (std::isnan(sample) || std::isinf(sample))
                        sample = 0.f;

                    samples[s] = sample;
                }
            }
        }

        // POST PROCESSING
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

        // INTERPOLATION NumTables
        // ...

        // EXPORT
        void makeWav(const juce::String& name)
        {
            const auto numSamples = buffer.getNumSamples();
            auto samples = buffer.getWritePointer(0);
            for (auto s = 0; s < numSamples; ++s)
                samples[s] = juce::jlimit(-1.f, 1.f, samples[s]);

            const auto path = juce::File::getSpecialLocation(
                juce::File::SpecialLocationType::userDesktopDirectory
            ).getFullPathName() + "\\WaveTableCreator\\";
            juce::File file(path);
            file.createDirectory();
            const auto fullName = name + ".wav";
            file = juce::File(path + "\\" + fullName);
            if (file.existsAsFile())
                file.deleteFile();
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
                writer->writeFromAudioSampleBuffer(buffer, 0, numSamples);
        }

    protected:
        const float tableSizeInv;
        const int tableSize;

        int numTables;
        Buffer buffer;
    };
}

int main (int, char*)
{
    using namespace wtc;

    const int tableSize = 2048;
    const int numTables = 256;

    TableFunc func = [numTables](float x, int tableIdx)
    {
        auto p = (float)tableIdx / numTables;

        auto a = std::rint(256.f * p);

        auto S = [](float x, float N, float M)
        {
            auto sum = 0.f;
			for (auto n = 1.f; n < N; ++n)
                sum += std::sin(M * n * x * pi) / (n * pi);

            return sum / M;
        };

        std::array<float, 8> fib =
        {
			1.f, 1.f, 2.f, 3.f, 5.f, 8.f, 13.f, 21.f
		};
		
		return S(x, a, fib[1]) * S(x, a, fib[5]) + S(x, a, fib[3]) * S(x, a, fib[2]);
    };

    Creator creator(tableSize, numTables, func);

    creator.dcOffsetEachTable();
    creator.normalizeEachTable();
    
    creator.makeWav("analogSawFibMulti");

    return 0;
}
