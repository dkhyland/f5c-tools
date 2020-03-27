#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <emmintrin.h>
#include <stdint.h>
#include <sys/time.h>

const int32_t N = 1000000;

int32_t * sequential(int32_t* a, int32_t* b, int32_t* c){
    for (int i=0; i < N; i++){
        a[i] = b[i]+c[i];
    }
    return a;
}

int32_t * parallel(int32_t* a, int32_t* b, int32_t* c){
    __m128i *av, *bv, *cv;
    // int32_t xx[N];
    av = (__m128i*)a; // assume 16-byte aligned
    bv = (__m128i*)b; // assume 16-byte aligned
    cv = (__m128i*)c; // assume 16-byte aligned
    for (int i = 0; i < N/4; i++){
        av[i] = _mm_add_epi32(bv[i], cv[i]);
        _mm_stream_si32(a,_mm_cvtsi128_si32(av[i]));
    }
    return a;
}

float * parallel_float(float* a, float* b, float* c){
    __m128 *av, *bv, *cv;
    // int32_t xx[N];
    av = (__m128*)a; // assume 16-byte aligned
    bv = (__m128*)b; // assume 16-byte aligned
    cv = (__m128*)c; // assume 16-byte aligned
    for (int i = 0; i < N/4; i++){
        av[i] = _mm_add_ps(bv[i], cv[i]);
        _mm_stream_ps(a,av[i]);
    }
    return a;
}

void test_runtime(){
    //Init
    int32_t * answer1, *answer2;
    int32_t * a;
    int32_t b[N];
    int32_t c[N];
    a = malloc(N * sizeof(int32_t));
    int i;
    for(i=0; i<N; ++i){
        b[i] = i;
        c[i] = i;
    }

    //Timing
    struct timeval start, stop;
    double cpu_time_used;

    //Run
    gettimeofday(&start,NULL);
    answer1 = sequential(a,b,c);
    gettimeofday(&stop,NULL);
    cpu_time_used = ((double) (stop.tv_usec - start.tv_usec));

    //Print
    printf("Sequential: %f ms\n",cpu_time_used);
    // for(i=0; i<N; ++i){
    //     printf("%d: %d, ",i,a[i]);
    // }

    //Run
    gettimeofday(&start,NULL);
    answer2 = parallel(a,b,c);
    gettimeofday(&stop,NULL);
    cpu_time_used = ((double) (stop.tv_usec - start.tv_usec));

    //Print
    printf("\nParallel: %f ms\n",cpu_time_used);
    // for(i=0; i<N; ++i){
    //     printf("%d: %d, ",i,a[i]);
    // }
    printf("\n");

    //test same values
    for(i=0; i<N; ++i){
        if(answer1[i] != answer2[i]){
            printf("%d not same. 1: %d, 2: %d\n",i,answer1[i],answer2[i]);
        }
    }

    free(a);
}

int32_t * compare(__m128 vec1, __m128 vec2){
    int32_t *ans;
    __m128i comp = _mm_add_epi32((__m128i)_mm_cmpeq_ps(vec1,vec2),_mm_set1_epi32(1));
    ans = (int32_t *)(&comp);
    return ans;
}

void test_compare(){
    int32_t *a = compare(_mm_set_ps(0,0.5,1.0,1.5),_mm_set_ps(0,1.0,1.0,1.0));
    for(int i=0; i<4; i++){
        printf("%d: %d\n",i,a[i]);
    }
}

void test_parallel_float(){
    //Init
    float * answer;
    float * a;
    float b[N];
    float c[N];
    a = malloc(N * sizeof(float));
    int i;
    for(i=0; i<N; ++i){
        b[i] = (float)i;
        c[i] = (float)i;
    }

    answer = parallel_float(a,b,c);
    for(i=0; i<N; ++i){
        printf("%d: %f\n",i,answer[i]);
    }
}

int main(){
    test_parallel_float();
}