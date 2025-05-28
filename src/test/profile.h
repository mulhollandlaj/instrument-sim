#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <math.h>

namespace profile {
    float lerp(float x0, float x1, float y0, float y1, float x) {
        return y0 + (x-x0)/(x1-x0) * (y1-y0);
    }

    template<int nx>
    struct Profile {
        std::array<float, nx> S;
        float length = 0;

        Profile(std::string fpath) {
            std::vector<float> dRaw, minRaw, maxRaw;
            std::array<float, nx> minGrid, maxGrid;
            unsigned int size = 0;
            float lengthmm = 0;

            std::ifstream file(fpath, std::ios::binary);
            if (file) {
                float data;
                file.read(reinterpret_cast<char*>(&size), sizeof(size));
                file.read(reinterpret_cast<char*>(&lengthmm), sizeof(lengthmm));
                for (int i = 0; i < size; i++) {
                    file.read(reinterpret_cast<char*>(&data), sizeof(data));
                    dRaw.push_back(data);
                }

                for (int i = 0; i < size; i++) {
                    file.read(reinterpret_cast<char*>(&data), sizeof(data));
                    minRaw.push_back(data);
                }

                for (int i = 0; i < size; i++) {
                    file.read(reinterpret_cast<char*>(&data), sizeof(data));
                    maxRaw.push_back(data);
                }
                file.close();
            } else {
                std::cerr << "Error opening file." << std::endl;
            }

            float h = lengthmm / float (nx);
            length = lengthmm * 1e-3;
            
            float x = 0;
            int i0 = 0;
            int i1 = 1;
            for (int i = 0; i < nx; i++) {
                x = h * ((float) i + 0.5);
                while (x > minRaw[i1]) {
                    i0++;
                    i1++;
                }
                minGrid[i] = lerp(minRaw[i0], minRaw[i1], dRaw[i0], dRaw[i1], x);
            }

            i0 = 0;
            i1 = 1;
            for (int i = 0; i < nx; i++) {
                x = h * ((float) i + 0.5);
                while (x >= maxRaw[i1]) {
                    i0++;
                    i1++;
                }
                maxGrid[i] = lerp(maxRaw[i0], maxRaw[i1], dRaw[i0], dRaw[i1], x);
            }

            for (int i = 0; i < nx; i++) {
                S[i] = M_PI*minGrid[i]*maxGrid[i]*1e-6*0.25;
            }
        }
    };
}