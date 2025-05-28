#include <fstream>
#include <algorithm>
#include <cstdint>
#include <string>

struct WavHeader {
    // Master RIFF chunk
    uint8_t const FileTypeBlocID[4] = {'R', 'I', 'F', 'F'};      // identifier 'RIFF'
    uint32_t mutable FileSize;                       // overall file size minus 8B
    uint8_t const FileFormat[4] = {'W', 'A', 'V', 'E'};          // identifier 'WAVE'

    // Chunk describing the data format
    uint8_t const FormatBlocID[4] = {'f', 'm', 't', ' '};        // identifier 'fmt '
    uint32_t const BlocSize = 16;                    // Chunk size minus 8B
    uint16_t const AudioFormat;                    // Audio format (1: PCM integer, 3: IEEE 754 float)
    uint16_t const NbrChannels;                    // Number of channels
    uint32_t const Frequency;                        // Sample rate (in Hertz)
    uint32_t const BytePerSec;                       // Number of bytes to read per second (Frequency * BytePerBloc)
    uint16_t const BytePerBloc;                    // Number of bytes per block (NbrChannels * BitsPerSample / 8)
    uint16_t const BitsPerSample;                  // Number of bits per sample

    // Chunk containing the sampled data
    uint8_t const DataBlocID[4] = {'d', 'a', 't', 'a'};          // identifier 'data'
    uint32_t mutable DataSize;                       // Sampled data size

    WavHeader (short AudioFormat, short NbrChannels, int Frequency, short BitsPerSample) : 
        FileSize(36),
        AudioFormat(AudioFormat),
        NbrChannels(NbrChannels),
        Frequency(Frequency),
        BitsPerSample(BitsPerSample),
        BytePerBloc(NbrChannels*BitsPerSample/8),
        DataSize(0),
        BytePerSec(Frequency*NbrChannels*BitsPerSample/8) {}
};


class WavWriter {
    private:
    std::ofstream stream;
    WavHeader* header;
    char const* outFile;

    template <typename T>
    void write(const T& t) {
        stream.write((const char*) &t, sizeof(T));
    }

    public:
    WavWriter(std::string outFile, short AudioFormat, short NbrChannels, int Frequency, short BitsPerSample) {
        stream = std::ofstream(outFile, std::ios::binary);
        header = new WavHeader(AudioFormat, NbrChannels, Frequency, BitsPerSample);
        write(header->FileTypeBlocID);
        write(header->FileSize);
        write(header->FileFormat);
        write(header->FormatBlocID);
        write(header->BlocSize);
        write(header->AudioFormat);
        write(header->NbrChannels);
        write(header->Frequency);
        write(header->BytePerSec);
        write(header->BytePerBloc);
        write(header->BitsPerSample);
        write(header->DataBlocID);
        write(header->DataSize);
    }
    
    void writeData(char* pBuffer, int size) {
        stream.write((char*) pBuffer, size);
        header->FileSize += size;
        header->DataSize += size;
    }

    void close() {
        stream.seekp(4);
        write(header->FileSize);
        stream.seekp(40);
        write(header->DataSize);
        stream.flush();
        stream.close();
    }
};