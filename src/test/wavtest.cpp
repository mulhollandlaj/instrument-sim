#include "wavwriter.h"
#include <math.h>
#include <string>

#define SAMPLE_RATE 44100
#define SIN_FREQ 415
#define AMPL 0x80000000

#define BUFFER_SIZE 512

#define TIME_SECONDS 5

int getSample(float t) {
    return (int) (AMPL*(
        0.001 * std::sin(M_PI*SIN_FREQ*t)
      + 0.01 * std::sin(5/4*M_PI*SIN_FREQ*t) 
      + 0.02 * std::sin(3/2*M_PI*SIN_FREQ*t)
      + 0.03 * std::sin(7/4*M_PI*SIN_FREQ*t)
      + 0.04 * std::sin(2*M_PI*SIN_FREQ*t)));
}

int main() {
    WavWriter ww = WavWriter((std::string) "test.wav", 1, 1, SAMPLE_RATE, 32);

    int sbuffer[BUFFER_SIZE];
    float t = 0;
    
    for (int i = 0; i < SAMPLE_RATE * TIME_SECONDS / BUFFER_SIZE; i++) {
        for (int j = 0; j < BUFFER_SIZE; j++) {
            sbuffer[j] = getSample(t);
            t += 1.0 / SAMPLE_RATE;
        }
        ww.writeData((char *) sbuffer, BUFFER_SIZE * sizeof(sbuffer[0]) / sizeof(char));
    }
    ww.close();
}