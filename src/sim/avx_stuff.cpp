#include <immintrin.h>
#include <algorithm>

template <size_t size>
class Avx512_floatarr {
    size_t blockCount = (size+15)/16;
    __m512 blockArr[blockCount];
    public:
        Avx512_floatarr(float * arr) {
            float arrCopy[blockCount*16];
            std::copy(arr, arr+size, arrCopy);
            for (size_t i = 0; i < blockCount; i++) {
                blockArr[i] = _mm512_set_ps(
                    arrCopy[16*i+15],
                    arrCopy[16*i+14],
                    arrCopy[16*i+13],
                    arrCopy[16*i+12],
                    arrCopy[16*i+11],
                    arrCopy[16*i+10],
                    arrCopy[16*i+9],
                    arrCopy[16*i+8],
                    arrCopy[16*i+7],
                    arrCopy[16*i+6],
                    arrCopy[16*i+5],
                    arrCopy[16*i+4],
                    arrCopy[16*i+3],
                    arrCopy[16*i+2],
                    arrCopy[16*i+1],
                    arrCopy[16*i]
                )
            }
        }
        
        float[size] toArray() {
            
        }
};