#include "f5c.h"
#include <emmintrin.h>
#include <assert.h>
#include <sys/time.h>
//#define DEBUG_ESTIMATED_SCALING 1
//#define DEBUG_RECALIB_SCALING 1
//#define DEBUG_ADAPTIVE 1

//Code was adapted from Nanopolish : nanopolish_raw_loader.cpp

//todo : can make more efficient using bit encoding
static inline uint32_t get_rank(char base) {
    if (base == 'A') { //todo: do we neeed simple alpha?
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else {
        WARNING("A None ACGT base found : %c", base);
        return 0;
    }
}

// return the lexicographic rank of the kmer amongst all strings of
// length k for this alphabet
static inline uint32_t get_kmer_rank(const char* str, uint32_t k) {
    //uint32_t p = 1;
    uint32_t r = 0;

    // from last base to first
    for (uint32_t i = 0; i < k; ++i) {
        //r += rank(str[k - i - 1]) * p;
        //p *= size();
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

//copy a kmer from a reference
static inline void kmer_cpy(char* dest, char* src, uint32_t k) {
    uint32_t i = 0;
    for (i = 0; i < k; i++) {
        dest[i] = src[i];
    }
    dest[i] = '\0';
}

scalings_t estimate_scalings_using_mom(char* sequence, int32_t sequence_len,
                                       model_t* pore_model, event_table et) {
    scalings_t out;
    int32_t n_kmers =
        sequence_len - KMER_SIZE + 1; //todo :strlen can be pre-calculated

    //const Alphabet* alphabet = pore_model.pmalphabet;

    // Calculate summary statistics over the events and
    // the model implied by the read
    float event_level_sum = 0.0f; //do we need double?
    for (uint32_t i = 0; i < et.n; ++i) {
        event_level_sum += et.event[i].mean;
    }

    float kmer_level_sum = 0.0f;
    float kmer_level_sq_sum = 0.0f;
    for (int32_t i = 0; i < n_kmers; ++i) {
        int32_t kr = get_kmer_rank(&sequence[i], KMER_SIZE);
        float l = pore_model[kr].level_mean;
        //fprintf(stderr,"Kmer : %c%c%c%c%c%c, kmer_rank : %d , kmer_mean : %f \n",sequence[i],sequence[i+1],sequence[i+2],sequence[i+3],sequence[i+4],sequence[i+5],kr,l);
        kmer_level_sum += l;
        kmer_level_sq_sum += l * l;
    }

    float shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

    // estimate scale
    float event_level_sq_sum = 0.0f;
    for (uint32_t i = 0; i < et.n; ++i) {
        event_level_sq_sum +=
            (et.event[i].mean - shift) * (et.event[i].mean - shift);
    }

    float scale = (event_level_sq_sum / et.n) / (kmer_level_sq_sum / n_kmers);

    //out.set4(shift, scale, 0.0, 1.0);
    out.shift = (float)shift;
    out.scale = (float)scale;

#ifdef DEBUG_ESTIMATED_SCALING
    fprintf(stderr, "event mean: %.2lf kmer mean: %.2lf shift: %.2lf\n",
            event_level_sum / et.n, kmer_level_sum / n_kmers, out.shift);
    fprintf(stderr, "event sq-mean: %.2lf kmer sq-mean: %.2lf scale: %.2lf\n",
            event_level_sq_sum / et.n, kmer_level_sq_sum / n_kmers, out.scale);
    //fprintf(stderr, "truth shift: %.2lf scale: %.2lf\n", pore_model.shift, pore_model.scale);
#endif
    return out;
}

static inline float log_normal_pdf(float x, float gp_mean, float gp_stdv,
                                   float gp_log_stdv) {
    /*INCOMPLETE*/
    float log_inv_sqrt_2pi = -0.918938f; // Natural logarithm
    float a = (x - gp_mean) / gp_stdv;
    return log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
    // return 1;
}

static inline float log_probability_match_r9(scalings_t scaling,
                                             model_t* models,
                                             event_table events, int32_t event_idx,
                                             int32_t kmer_rank, uint8_t strand,
                                             float sample_rate) {
    // event level mean, scaled with the drift value
    strand = 0;
    assert(kmer_rank < 4096);
    //float level = read.get_drift_scaled_level(event_idx, strand);

    //float time =
    //    (events.event[event_idx].start - events.event[0].start) / sample_rate;
    float unscaledLevel = events.event[event_idx].mean;
    float scaledLevel = unscaledLevel;
    //float scaledLevel = unscaledLevel - time * scaling.shift;

    //fprintf(stderr, "level %f\n",scaledLevel);
    //GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    float gp_mean =
        scaling.scale * models[kmer_rank].level_mean + scaling.shift;
    float gp_stdv = models[kmer_rank].level_stdv * 1; //scaling.var = 1;
    // float gp_stdv = 0;
    // float gp_log_stdv = models[kmer_rank].level_log_stdv + scaling.log_var;
    // if(models[kmer_rank].level_stdv <0.01 ){
    // 	fprintf(stderr,"very small std dev %f\n",models[kmer_rank].level_stdv);
    // }
    #ifdef CACHED_LOG
        float gp_log_stdv = models[kmer_rank].level_log_stdv;
    #else
        float gp_log_stdv =
        log(models[kmer_rank].level_stdv); // scaling.log_var = log(1)=0;
    #endif

    float lp = log_normal_pdf(scaledLevel, gp_mean, gp_stdv, gp_log_stdv);
    return lp;
}

//SSE2 helper function to compare two floating point vectors. Returns 1 in positions where the value is the same and 0 if not.
__m128i compare_from_vec(__m128 vec1, __m128 vec2){
    return _mm_add_epi32((__m128i)_mm_cmpeq_ps(vec1,vec2),_mm_set1_epi32(1));
}

//Helper function to print a single band vector and a from vector
void print_band_trace(int32_t band_idx, __m128 band_vec, __m128i from_vec){
    float * band = (float *)malloc(4 * sizeof(float));
    int32_t * from = (int32_t *)malloc(4 * sizeof(int32_t));;

    _mm_store_ps(band,band_vec);
    //_mm_stream_si32(from,_mm_cvtsi128_si32(from_vec));
    from[0] = _mm_cvtsi128_si32(from_vec);
    from[1] = _mm_cvtsi128_si32(_mm_srai_epi32(from_vec,32));
    from[2] = _mm_cvtsi128_si32(_mm_srai_epi32(from_vec,64));
    from[3] = _mm_cvtsi128_si32(_mm_srai_epi32(from_vec,96));

    // fprintf(stderr, "Band number: %d, Band: (%.2f,%.2f,%.2f,%.2f), From: (%d,%d,%d,%d)\n",band_idx,band[0],band[1],band[2],band[3],from[0],from[1],from[2],from[3]);


    free(band);
    free(from);
}

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)

#define move_down(curr_band) \
    { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) \
    { curr_band.event_idx, curr_band.kmer_idx + 1 }

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


#ifdef ALIGN_2D_ARRAY
    #define BAND_ARRAY(r, c) ( bands[(r)][(c)] )
    #define TRACE_ARRAY(r, c) ( trace[(r)][(c)] )
#else
    #define BAND_ARRAY(r, c) ( bands[((r)*(ALN_BANDWIDTH)+(c))] )
    #define TRACE_ARRAY(r, c) ( trace[((r)*(ALN_BANDWIDTH)+(c))] )
#endif

//SIMD INSTRUCTIONS. TODO: implement

//Return the vector corresponding to the (r,c) position in band/trace vec
#define BAND_ARRAY_VEC(r,c) ( band_vec[(r)*(ALN_BANDWIDTH/4)+(c)] )
#define TRACE_ARRAY_VEC(r,c) ( trace_vec[(r)*(ALN_BANDWIDTH/4)+(c)] )
#define sse2_convert_size(size) ((size) + 3 / 4)
#define sse2_convert_index(index) ((index) / 4)

int32_t align_simd(AlignedPair* out_2, char* sequence, int32_t sequence_len,
              event_table events, model_t* models, scalings_t scaling,
              float sample_rate) {
    //fprintf(stderr, "%s\n", sequence);
    //fprintf(stderr, "Scaling %f %f", scaling.scale, scaling.shift);

    int32_t strand_idx = 0;
    int32_t k = 6;

    // int32_t n_events = events[strand_idx].n;
    int32_t n_events = events.n;
    int32_t n_kmers = sequence_len - k + 1;
    //fprintf(stderr,"n_kmers : %d\n",n_kmers);
    // backtrack markers
    const int32_t FROM_D = 0;
    const int32_t FROM_U = 1;
    const int32_t FROM_L = 2;

    // qc
    float min_average_log_emission = -5.0;
    int max_gap_threshold = 50;

    // banding
    int32_t bandwidth = ALN_BANDWIDTH;
    int32_t half_bandwidth = ALN_BANDWIDTH / 2;

    //SSE2 number of vectors in a band
    // int32_t bandwidth_vec = sse2_convert_size(bandwidth);

    // transition penalties
    float events_per_kmer = (float)n_events / n_kmers;
    float p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    float epsilon = 1e-10;
    float lp_skip = log(epsilon);
    float lp_stay = log(p_stay);
    float lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    float lp_trim = log(0.01);

    //SSE2 constant vectors
    __m128 lp_skip_vec = _mm_set1_ps(lp_skip);
    __m128 lp_stay_vec = _mm_set1_ps(lp_stay);
    __m128 lp_step_vec = _mm_set1_ps(lp_step);

    // dp matrix
    int32_t n_rows = n_events + 1;
    int32_t n_cols = n_kmers + 1;
    int32_t n_bands = n_rows + n_cols;

    // Initialize

    // Precompute k-mer ranks to avoid doing this in the inner loop
    int32_t* kmer_ranks = (int32_t*)malloc(sizeof(int32_t) * n_kmers);
    MALLOC_CHK(kmer_ranks);

    for (int32_t i = 0; i < n_kmers; ++i) {
        //>>>   >>>>>> New replacement begin
        char* substring = &sequence[i];
        kmer_ranks[i] = get_kmer_rank(substring, k);
        //<<<<<<<<< New replacement over
    }

#ifdef ALIGN_2D_ARRAY
    float** bands = (float**)malloc(sizeof(float*) * n_bands);
    MALLOC_CHK(bands);
    int32_t** trace = (int32_t**)malloc(sizeof(int32_t*) * n_bands);
    MALLOC_CHK(trace);
#else
    float* bands = (float*)malloc(sizeof(float) * n_bands * bandwidth);
    MALLOC_CHK(bands);
    int32_t* trace = (int32_t*)malloc(sizeof(int32_t) * n_bands * bandwidth);
    MALLOC_CHK(trace);

    //SSE2 version of bands and trace
    // __m128 *band_vec = (__m128 *)malloc(sizeof(__m128) * n_bands * bandwidth_vec);
    // MALLOC_CHK(band_vec);
    // __m128i *trace_vec = (__m128i *)malloc(sizeof(__m128i) * n_bands * bandwidth_vec);
    // MALLOC_CHK(band_vec);

#endif

    //Initialize default values
    for (int32_t i = 0; i < n_bands; i++) {
    #ifdef ALIGN_2D_ARRAY
        bands[i] = (float*)malloc(sizeof(float) * bandwidth);
        MALLOC_CHK(bands[i]);
        trace[i] = (int32_t*)malloc(sizeof(int32_t) * bandwidth);
        MALLOC_CHK(trace[i]);
    #endif

        for (int32_t j = 0; j < bandwidth; j++) {
            BAND_ARRAY(i,j) = -INFINITY;
            TRACE_ARRAY(i,j) = 0;
        }
    }

    //SSE2 initialize default values
    // for (int32_t i = 0; i < n_bands; i++) {
    //     for (int32_t j = 0; j < bandwidth_vec; j++) {
    //         BAND_ARRAY_VEC(i,j) = _mm_set1_ps(-INFINITY);
    //         TRACE_ARRAY_VEC(i,j) = _mm_set1_epi32(0);
    //     }
    // }

    // Keep track of the event/kmer index for the lower left corner of the band
    // these indices are updated at every iteration to perform the adaptive banding
    // Only the first two  have their coordinates initialized, the rest are computed adaptively

    struct EventKmerPair {
        int32_t event_idx;
        int32_t kmer_idx;
    };
    //>>>>>>>>>>>>>>>>>New Replacement Begin
    struct EventKmerPair* band_lower_left =
        (struct EventKmerPair*)malloc(sizeof(struct EventKmerPair) * n_bands);
    MALLOC_CHK(band_lower_left);
    //std::vector<EventKmerPair> band_lower_left(n_);
    //<<<<<<<<<<<<<<<<<New Replacement over

    // initialize range of first two
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    //simd_debug
    // fprintf(stderr,"init: %d,%d,%d,%d\n",band_lower_left[0].event_idx,band_lower_left[0].kmer_idx,
    // band_lower_left[1].event_idx,band_lower_left[1].kmer_idx);

    int32_t start_cell_offset = band_kmer_to_offset(0, -1);
    // assert(is_offset_valid(start_cell_offset));
    // assert(band_event_to_offset(0, -1) == start_cell_offset);
    BAND_ARRAY(0,start_cell_offset) = 0.0f;

    //Put the zero in the right place within the vector
    // int32_t sse2_start_position = start_cell_offset % 4;
    // int32_t sse2_start_cell = sse2_convert_index(start_cell_offset);
    // switch (sse2_start_position)
    // {
    //     case 0:
    //         BAND_ARRAY_VEC(0,sse2_start_cell) = _mm_set_ps(0.0,-INFINITY,-INFINITY,-INFINITY);
    //     case 1:
    //         BAND_ARRAY_VEC(0,sse2_start_cell) = _mm_set_ps(-INFINITY,0.0,-INFINITY,-INFINITY);
    //     case 2:
    //         BAND_ARRAY_VEC(0,sse2_start_cell) = _mm_set_ps(-INFINITY,-INFINITY,0.0,-INFINITY);
    //     case 3:
    //         BAND_ARRAY_VEC(0,sse2_start_cell) = _mm_set_ps(-INFINITY,-INFINITY,-INFINITY,0.0);
    // }

    // band 1: first event is trimmed
    int32_t first_trim_offset = band_event_to_offset(1, 0);
    // assert(kmer_at_offset(1, first_trim_offset) == -1);
    // assert(is_offset_valid(first_trim_offset));
    BAND_ARRAY(1,first_trim_offset) = lp_trim;
    TRACE_ARRAY(1,first_trim_offset) = FROM_U;

    //Put the values in the right place within the vector
    // int32_t sse2_first_trim = sse2_convert_index(first_trim_offset);
    // switch (sse2_start_position)
    // {
    //     case 0:
    //         BAND_ARRAY_VEC(1,sse2_first_trim) = _mm_set_ps(lp_trim,-INFINITY,-INFINITY,-INFINITY);
    //         TRACE_ARRAY_VEC(1,sse2_first_trim) = _mm_set_epi32(FROM_U,0,0,0);
    //     case 1:
    //         BAND_ARRAY_VEC(1,sse2_first_trim) = _mm_set_ps(-INFINITY,lp_trim,-INFINITY,-INFINITY);
    //         TRACE_ARRAY_VEC(1,sse2_first_trim) = _mm_set_epi32(0,FROM_U,0,0);
    //     case 2:
    //         BAND_ARRAY_VEC(1,sse2_first_trim) = _mm_set_ps(-INFINITY,-INFINITY,lp_trim,-INFINITY);
    //         TRACE_ARRAY_VEC(1,sse2_first_trim) = _mm_set_epi32(0,0,FROM_U,0);
    //     case 3:
    //         BAND_ARRAY_VEC(1,sse2_first_trim) = _mm_set_ps(-INFINITY,-INFINITY,-INFINITY,lp_trim);
    //         TRACE_ARRAY_VEC(1,sse2_first_trim) = _mm_set_epi32(0,0,0,FROM_U);
    // }

    //Declare/initialise variables needed for the loops
    int32_t event_idx,kmer_idx,kmer_rank,offset_up,offset_left,offset_diag;
    float lp_emission;
    int32_t fills = 0;
    float *up_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(up_arr);
    float *left_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(left_arr);
    float *diag_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(diag_arr);
    float *lp_emission_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(lp_emission_arr);
    int32_t *from_arr = (int32_t *)malloc(sizeof(int32_t) * 4);
    MALLOC_CHK(from_arr);
    float * band_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(band_arr);

#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1,
            first_trim_offset, 0, -1, bands[1][first_trim_offset]);
#endif
    //simd_debug
    int num_right = 0;
    int num_down = 0;
    // fill in remaining bands
    for (int32_t band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
        float ll = BAND_ARRAY(band_idx - 1,0);
        float ur = BAND_ARRAY(band_idx - 1,bandwidth - 1);

        //simd_debug
        fprintf(stderr,"band_idx: %d, ll: %f,ur: %f\n",band_idx,ll,ur);

        bool ll_ob = ll == -INFINITY;
        bool ur_ob = ur == -INFINITY;

        bool right = false;
        if (ll_ob && ur_ob) {
            right = band_idx % 2 == 1;
        } else {
            right = ll < ur; // Suzuki's rule
        }

        if (right) {
            band_lower_left[band_idx] =
                move_right(band_lower_left[band_idx - 1]);
            //simd_debug
            num_right++;
        } else {
            band_lower_left[band_idx] =
                move_down(band_lower_left[band_idx - 1]);
            //simd_debug
            num_down++;
        }
        // If the trim state is within the band, fill it in here
        int32_t trim_offset = band_kmer_to_offset(band_idx, -1);
        if (is_offset_valid(trim_offset)) {
            int32_t event_idx = event_at_offset(band_idx, trim_offset);
            if (event_idx >= 0 && event_idx < n_events) {
                BAND_ARRAY(band_idx,trim_offset) = lp_trim * (event_idx + 1);
                TRACE_ARRAY(band_idx,trim_offset) = FROM_U;
            } else {
                BAND_ARRAY(band_idx,trim_offset) = -INFINITY;
            }
        }

        // Get the offsets for the first and last event and kmer
        // We restrict the inner loop to only these values
        int32_t kmer_min_offset = band_kmer_to_offset(band_idx, 0);
        int32_t kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
        int32_t event_min_offset = band_event_to_offset(band_idx, n_events - 1);
        int32_t event_max_offset = band_event_to_offset(band_idx, -1);

        int32_t min_offset = MAX(kmer_min_offset, event_min_offset);
        min_offset = MAX(min_offset, 0);

        int32_t max_offset = MIN(kmer_max_offset, event_max_offset);
        max_offset = MIN(max_offset, bandwidth);

        //Inner loop: Parallelised with SSE2. Jump up 4 every time
        for (int32_t offset = min_offset; offset < max_offset; offset += 4) {
            if(offset + 4 >= max_offset){
                //If we don't have >= 4 cells left to fill, compute individually using a sequential loop
                for(int32_t seq_offset = offset; seq_offset < max_offset; seq_offset++){
                    event_idx = event_at_offset(band_idx, seq_offset);
                    kmer_idx = kmer_at_offset(band_idx, seq_offset);
                    kmer_rank = kmer_ranks[kmer_idx];
                    //simd_debug
                    // fprintf(stderr, "event idx %d, kmer_idx: %d, kmer rank %d, band_idx: %d, seq_offset: %d, bllevent: %d, bllkmer: %d\n", 
                    // event_idx,kmer_idx,kmer_rank,band_idx,seq_offset,band_lower_left[band_idx].event_idx,band_lower_left[band_idx].kmer_idx);

                    //Offset of the up, left, and diagonal positions
                    offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
                    offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
                    offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
                    // verify loop conditions
                    assert(kmer_idx >= 0 && kmer_idx < n_kmers);
                    assert(event_idx >= 0 && event_idx < n_events);
                    assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
                    assert(offset_up - offset_left == 1);
                    assert(seq_offset >= 0 && seq_offset < bandwidth);
#endif

                    float up = is_offset_valid(offset_up)
                            ? BAND_ARRAY(band_idx - 1,offset_up)
                            : -INFINITY;

                    float left = is_offset_valid(offset_left)
                                ? BAND_ARRAY(band_idx - 1,offset_left)
                                : -INFINITY;

                    float diag = is_offset_valid(offset_diag)
                                ? BAND_ARRAY(band_idx - 2,offset_diag)
                                : -INFINITY;

                    lp_emission = log_probability_match_r9(scaling, models, events, event_idx,
                                            kmer_rank, strand_idx, sample_rate);
                    
                    //simd_debug
                    // fprintf(stderr, "lp emission : %f\n", lp_emission);
                    
                    //Compute score and from values for single entry
                    float score_d_s = diag + lp_step + lp_emission;
                    float score_u_s = up + lp_stay + lp_emission;
                    float score_l_s = left + lp_skip;

                    //A single max_score/from calculation
                    float max_score_single = score_d_s;
                    int32_t from_single = FROM_D;

                    max_score_single = score_u_s > max_score_single ? score_u_s : max_score_single;
                    from_single = max_score_single == score_u_s ? FROM_U : from_single;
                    max_score_single = score_l_s > max_score_single ? score_l_s : max_score_single;
                    from_single = max_score_single == score_l_s ? FROM_L : from_single;

                    //simd_debug
                    fprintf(stderr, "offset: %d, score: %f, from : %d\n",seq_offset,max_score_single,from_single);

                    //Store in arrays
                    BAND_ARRAY(band_idx,seq_offset) = max_score_single;
                    TRACE_ARRAY(band_idx,seq_offset) = from_single;
                    fills += 1;
                }
            }else{
                //Compute using SIMD
                //Need to load values sequentially because the __m128i vectors don't overlap
                //Load 4 values corresponding to the left, up and diagonal bands into arrays
                for (int32_t vec_pos = 0; vec_pos < 4 ; ++vec_pos) {
                    //Index of the first element of the vector we are looking at
                    event_idx = event_at_offset(band_idx, offset + vec_pos);
                    kmer_idx = kmer_at_offset(band_idx, offset + vec_pos);
                    //simd_debug
                    // fprintf(stderr, "event idx %d, kmer_idx: %d, kmer rank %d, band_idx: %d, vec_pos: %d, offset: %d, bllevent: %d, bllkmer: %d\n", event_idx,kmer_idx,kmer_rank,
                    // band_idx,vec_pos,offset,band_lower_left[band_idx].event_idx,band_lower_left[band_idx].kmer_idx);
                    kmer_rank = kmer_ranks[kmer_idx];

                    //Offset of the up, left, and diagonal positions
                    offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
                    offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
                    offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
                    // verify loop conditions
                    assert(kmer_idx >= 0 && kmer_idx < n_kmers);
                    assert(event_idx >= 0 && event_idx < n_events);
                    assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
                    assert(offset_up - offset_left == 1);
                    assert(offset >= 0 && offset < bandwidth);
#endif

                    float up = is_offset_valid(offset_up)
                            ? BAND_ARRAY(band_idx - 1,offset_up)
                            : -INFINITY;
                    up_arr[vec_pos] = up;

                    float left = is_offset_valid(offset_left)
                                ? BAND_ARRAY(band_idx - 1,offset_left)
                                : -INFINITY;
                    left_arr[vec_pos] = left;

                    float diag = is_offset_valid(offset_diag)
                                ? BAND_ARRAY(band_idx - 2,offset_diag)
                                : -INFINITY;
                    diag_arr[vec_pos] = diag;

                    lp_emission = log_probability_match_r9(scaling, models, events, event_idx,
                                            kmer_rank, strand_idx, sample_rate);
                    lp_emission_arr[vec_pos] = lp_emission;
                    //simd_debug
                    //fprintf(stderr, "lp emission : %f\n", lp_emission);
                }

                //convert data from the arrays to __m128
                __m128 up_vec = _mm_set_ps(up_arr[0],up_arr[1],up_arr[2],up_arr[3]);
                __m128 left_vec = _mm_set_ps(left_arr[0],left_arr[1],left_arr[2],left_arr[3]);
                __m128 diag_vec = _mm_set_ps(diag_arr[0],diag_arr[1],diag_arr[2],diag_arr[3]);
                __m128 lp_emission_vec = _mm_set_ps(lp_emission_arr[0],lp_emission_arr[1],lp_emission_arr[2],lp_emission_arr[3]);

                __m128 score_d = _mm_add_ps(diag_vec,_mm_add_ps(lp_step_vec,lp_emission_vec));
                __m128 score_u = _mm_add_ps(up_vec,_mm_add_ps(lp_stay_vec,lp_emission_vec));
                __m128 score_l = _mm_add_ps(left_vec,lp_skip_vec);

                __m128 max_score = score_d;
                max_score = _mm_max_ps(score_l,_mm_max_ps(score_u,max_score));

                //These vectors have a 1 where the max_score corresponds to the direction, and 0 where it doesn't
                __m128i compare_up = compare_from_vec(max_score,score_u);
                __m128i compare_left = compare_from_vec(max_score,score_l);
                //FROM_D=0, FROM_U=1, FROM_L=2, so only need to add compare_up to 2 * compare_left
                __m128i from = _mm_add_epi32(compare_up,_mm_add_epi32(compare_left,compare_left));

                // print_band_trace(band_idx,max_score,from);

                //Store results in BAND and TRACE array
                float *band_position = &(BAND_ARRAY(band_idx,offset));
                int32_t *trace_position = &(TRACE_ARRAY(band_idx,offset));

                _mm_storer_ps(band_arr,max_score); 
                band_position[0] = band_arr[0];
                band_position[1] = band_arr[1];
                band_position[2] = band_arr[2];
                band_position[3] = band_arr[3];
                
                trace_position[0] = _mm_cvtsi128_si32(from);     
                trace_position[1] = _mm_cvtsi128_si32(_mm_srli_epi32(from,32));
                trace_position[2] = _mm_cvtsi128_si32(_mm_srli_epi32(from,64));
                trace_position[3] = _mm_cvtsi128_si32 (_mm_srli_epi32(from,96));
                

                //simd_debug
                //fprintf(stderr,"bands: %f %f %f %f\n",band_position[0],band_position[1],band_position[2],band_position[3]);
                fprintf(stderr,"band array: %f %f %f %f\n",BAND_ARRAY(band_idx,offset),BAND_ARRAY(band_idx,offset+1),BAND_ARRAY(band_idx,offset+2),BAND_ARRAY(band_idx,offset+3));
                //simd_debug
                //fprintf(stderr,"traces: %d %d %d %d\n",trace_position[0],trace_position[1],trace_position[2],trace_position[3]);
                fprintf(stderr,"offset: %d, trace array: %d %d %d %d\n",offset,TRACE_ARRAY(band_idx,offset),TRACE_ARRAY(band_idx,offset+1),TRACE_ARRAY(band_idx,offset+2),TRACE_ARRAY(band_idx,offset+3));
                
                fills += 4;
            }
#ifdef DEBUG_ADAPTIVE
            fprintf(stderr,
                    "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n",
                    offset_up, offset_diag, offset_left);
            fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up,
                    diag, left);
            fprintf(stderr,
                    "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: "
                    "%.2lf\n",
                    band_idx, offset, event_idx, kmer_idx, max_score, from,
                    lp_emission);
#endif
        }
    }
    //simd_debug
    // fprintf(stderr,"right: %d, down: %d\n",num_right,num_down);
    //
    // Backtrack to compute alignment
    //
    double sum_emission = 0;
    double n_aligned_events = 0;

    //>>>>>>>>>>>>>> New replacement begin
    // std::vector<AlignedPair> out;

    int32_t outIndex = 0;
    //<<<<<<<<<<<<<<<<New Replacement over

    float max_score = -INFINITY;
    int32_t curr_event_idx = 0;
    int32_t curr_kmer_idx = n_kmers - 1;

    // Find best score between an event and the last k-mer. after trimming the remaining evnets
    for (int32_t event_idx = 0; event_idx < n_events; ++event_idx) {
        int32_t band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);

        //>>>>>>>New  replacement begin
        /*assert(band_idx < bands.size());*/
        assert((int32_t)band_idx < n_bands);
        //<<<<<<<<New Replacement over
        int32_t offset = band_event_to_offset(band_idx, event_idx);

        //simd_debug
        if(event_idx == 0){
        fprintf(stderr,"bll[bi]: %d, offset: %d, band_idx: %d, event_idx: %d\n",band_lower_left[band_idx].event_idx,offset,band_idx,event_idx);
        }

        if (is_offset_valid(offset)) {
            float s =
                BAND_ARRAY(band_idx,offset) + (n_events - event_idx) * lp_trim;
            //simd_debug
            // fprintf(stderr,"s: %f, offset: %d, band_idx: %d, event_idx: %d\n",s,offset,band_idx,event_idx);
            if (s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }
    //simd_debug
    // fprintf(stderr,"max_score: %f, n_events: %d, curr_kmer_idx: %d\n",max_score,n_events,curr_kmer_idx);

#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx,
            curr_kmer_idx, max_score);
#endif

    int32_t curr_gap = 0;
    int32_t max_gap = 0;
    while (curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        // emit alignment
        //>>>>>>>New Repalcement begin
        assert(outIndex < (int32_t)(n_events * 2));
        out_2[outIndex].ref_pos = curr_kmer_idx;
        out_2[outIndex].read_pos = curr_event_idx;
        outIndex++;
        // out.push_back({curr_kmer_idx, curr_event_idx});
        //<<<<<<<<<New Replacement over

#ifdef DEBUG_ADAPTIVE
        fprintf(stderr, "[adaback] ei: %d ki: %d\n", curr_event_idx,
                curr_kmer_idx);
#endif
        // qc stats
        //>>>>>>>>>>>>>>New Replacement begin
        char* substring = &sequence[curr_kmer_idx];
        int32_t kmer_rank = get_kmer_rank(substring, k);
        //<<<<<<<<<<<<<New Replacement over
        float tempLogProb = log_probability_match_r9(
            scaling, models, events, curr_event_idx, kmer_rank, 0, sample_rate);

        //simd_debug
        // fprintf(stderr,"curr_event: %d, kmer_rank: %d, sample: %f\n",curr_event_idx,kmer_rank,sample_rate);

        sum_emission += tempLogProb;
        // fprintf(stderr, "lp_emission %f \n", tempLogProb);
        //simd_debug
        // fprintf(stderr,"lp_emission %f, sum_emission %f, n_aligned_events %d\n",tempLogProb,sum_emission,outIndex);

        n_aligned_events += 1;

        int32_t band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        int32_t offset = band_event_to_offset(band_idx, curr_event_idx);
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

        int32_t from = TRACE_ARRAY(band_idx,offset);
        if (from == FROM_D) {
            curr_kmer_idx -= 1;
            curr_event_idx -= 1;
            curr_gap = 0;
        } else if (from == FROM_U) {
            curr_event_idx -= 1;
            curr_gap = 0;
        } else {
            curr_kmer_idx -= 1;
            curr_gap += 1;
            max_gap = MAX(curr_gap, max_gap);
        }
    }

    //>>>>>>>>New replacement begin
    // std::reverse(out.begin(), out.end());
    int32_t c;
    int32_t end = outIndex - 1;
    for (c = 0; c < outIndex / 2; c++) {
        int32_t ref_pos_temp = out_2[c].ref_pos;
        int32_t read_pos_temp = out_2[c].read_pos;
        out_2[c].ref_pos = out_2[end].ref_pos;
        out_2[c].read_pos = out_2[end].read_pos;
        out_2[end].ref_pos = ref_pos_temp;
        out_2[end].read_pos = read_pos_temp;
        end--;
    }

    // if(outIndex>1){
    //   AlignedPair temp={out_2[0].ref_pos,out[0].read_pos};
    //   int i;
    //   for(i=0;i<outIndex-1;i++){
    //     out_2[i]={out_2[outIndex-1-i].ref_pos,out[outIndex-1-i].read_pos};
    //   }
    //   out[outIndex-1]={temp.ref_pos,temp.read_pos};
    // }
    //<<<<<<<<<New replacement over

    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
    //simd_debug
    // fprintf(stderr,"sum_emission %f, n_aligned_events %f, avg_log_emission %f\n",sum_emission,n_aligned_events,avg_log_emission);
    //>>>>>>>>>>>>>New replacement begin
    bool spanned = out_2[0].ref_pos == 0 &&
                   out_2[outIndex - 1].ref_pos == int32_t(n_kmers - 1);
    // bool spanned = out.front().ref_pos == 0 && out.back().ref_pos == n_kmers - 1;
    //<<<<<<<<<<<<<New replacement over
    //bool failed = false;

    //simd_debug
    if(avg_log_emission < min_average_log_emission){
        fprintf(stderr,"avg_log_emission < min_average_log_emission: %f, %f\n",avg_log_emission,min_average_log_emission);
    }
    //simd_debug
    if( !spanned ){
        fprintf(stderr,"not spanned\n");
    }
    //simd_debug
    if(max_gap > max_gap_threshold){
        fprintf(stderr,"max_gap > max_gap_threshold: %d, %d\n",max_gap,max_gap_threshold);
    }

    if (avg_log_emission < min_average_log_emission || !spanned ||
        max_gap > max_gap_threshold) {
        //failed = true;
        //>>>>>>>>>>>>>New replacement begin
        outIndex = 0;
        // out.clear();
        //free(out_2);
        //out_2 = NULL;
        //<<<<<<<<<<<<<New replacement over
    }

    free(kmer_ranks);
#ifdef ALIGN_2D_ARRAY
    for (int32_t i = 0; i < n_bands; i++) {
        free(bands[i]);
        free(trace[i]);
    }
#endif

    //free regular mallocs
    free(bands);
    free(trace);
    free(band_lower_left);

    //free sse2 mallocs
    // free(band_vec);
    // free(trace_vec);
    free(up_arr);
    free(left_arr);
    free(diag_arr);
    free(lp_emission_arr);
    free(from_arr);
    free(band_arr);

    //fprintf(stderr, "ada\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, max_gap, fills);
    //outSize=outIndex;
    //if(outIndex>500000)fprintf(stderr, "Max outSize %d\n", outIndex);
    return outIndex;
}

int32_t postalign(event_alignment_t* alignment, index_pair_t* base_to_event_map,double* events_per_base,
                  char* sequence, int32_t n_kmers, AlignedPair* event_alignment,
                  int32_t n_events) {
    /* transform alignment into the base-to-event map*/
    // create base-to-event map
    // index_pair_t* base_to_event_map =
    //     (index_pair_t*)(malloc(sizeof(index_pair_t) * n_kmers));
    // MALLOC_CHK(base_to_event_map);

    //initialisesing (todo : check if really required)
    int32_t i = 0;
    for (i = 0; i < n_kmers; i++) {
        base_to_event_map[i].start = -1;
        base_to_event_map[i].stop = -1;
    }

    int32_t max_event = 0;
    int32_t min_event = INT32_MAX;

    int32_t prev_event_idx = -1;

    for (i = 0; i < n_events; ++i) {
        int32_t k_idx = event_alignment[i].ref_pos;
        int32_t event_idx = event_alignment[i].read_pos;
        index_pair_t* elem = &base_to_event_map[k_idx];
        //fprintf(stderr, "eventpar %d %d k_idx %d event_idx %d\n",elem.start, elem.stop,k_idx,event_idx);
        if (event_idx != prev_event_idx) {
            if (elem->start == -1) {
                elem->start = event_idx;
            }
            elem->stop = event_idx;
        }
        max_event = max_event > event_idx ? max_event : event_idx;
        min_event = min_event < event_idx ? min_event : event_idx;
        prev_event_idx = event_idx;
    }

    // for (i = 0; i < n_kmers; ++i) {
    //     fprintf(stderr,"base_to_event_map - start %d stop %d\n", base_to_event_map[i].start,base_to_event_map[i].stop);
    // }

    *events_per_base = (float)(max_event - min_event) / n_kmers;

    /*prepare data structures for the final calibration*/

    int32_t alignment_index = 0;
    int32_t prev_kmer_rank = -1;

    int32_t ki;
    for (ki = 0; ki < n_kmers; ++ki) {
        index_pair_t event_range_for_kmer = base_to_event_map[ki];

        //fprintf(stderr, "kindex %d base_to_event_map - start %d stop %d\n",ki,event_range_for_kmer.start, event_range_for_kmer.stop);

        // skip kmers without events
        if (event_range_for_kmer.start == -1) {
            continue;
        }

        //skip k-mers that cannot be shifted to a valid position
        // int32_t shift_offset=0;
        // if(ki + shift_offset < 0 || ki + shift_offset >= n_kmers) {
        //     continue;
        // }

        for (int32_t event_idx = event_range_for_kmer.start;
             event_idx <= event_range_for_kmer.stop; event_idx++) {
            //fprintf(stderr,"event idx %d n events %d\n",event_idx,n_events);
            // assert(event_idx < n_events);

            // since we use the 1D read seqence here we never have to reverse complement
            int32_t kmer_rank = get_kmer_rank(&sequence[ki], KMER_SIZE);

            event_alignment_t ea;
            // ref data
            //ea.ref_name = "read";
            ea.read_idx = -1; // not needed
            kmer_cpy(ea.ref_kmer, &sequence[ki], KMER_SIZE);
            ea.ref_position = ki;
            //ea.strand_idx = strand_idx;
            ea.event_idx = event_idx;
            ea.rc = false;
            kmer_cpy(ea.model_kmer, &sequence[ki], KMER_SIZE);
            ea.hmm_state = prev_kmer_rank != kmer_rank ? 'M' : 'E';
            if (alignment_index >
                n_events) { //todo : this is ugly. check and fix.
                ERROR("We have run out of space in event_alignment_t* "
                      "alignment. Assumption fialed. Current size %d",
                      n_events);
                exit(EXIT_FAILURE);
            }
            alignment[alignment_index] = ea;
            alignment_index++;
            prev_kmer_rank = kmer_rank;
        }
        //fprintf(stderr,"event idx : %d\n",alignment_index);
    }

    //free(base_to_event_map);
    return alignment_index;
}

// recalculate shift, scale, drift, scale_sd from an alignment and the read
// returns true if the recalibration was performed
// in either case, sets residual to the L1 norm of the residual
bool recalibrate_model(model_t* pore_model, event_table et,
                       scalings_t* scallings,
                       const event_alignment_t* alignment_output,
                       int32_t num_alignments, bool scale_var) {
    //std::vector<double> raw_events, times, level_means, level_stdvs;
    //std::cout << "Previous pore model parameters: " << sr.pore_model[strand_idx].shift << ", "
    //                                                << sr.pore_model[strand_idx].scale << ", "
    //                                                << sr.pore_model[strand_idx].drift << ", "
    //                                                << sr.pore_model[strand_idx].var   << std::endl;

    // extract necessary vectors from the read and the pore model; note do not want scaled values
    int32_t num_M_state = 0;
    for (int32_t ei = 0; ei < num_alignments; ++ei) {
        event_alignment_t ea = alignment_output[ei];
        if (ea.hmm_state == 'M') {
            num_M_state++;
            //
            //fprintf(stdout, "recalibrate ei: %zu level: %.2lf kmer: %s model: %.2lf\n",
            //        ei, sr.get_uncorrected_level(ea.event_idx, strand_idx), model_kmer.c_str(),
            //        sr.pore_model[strand_idx].states[rank].level_mean);
            //
        }
    }

    const int32_t minNumEventsToRescale = 200;
    bool recalibrated = false;
    if (num_M_state >= minNumEventsToRescale) {
        // Assemble linear system corresponding to weighted least squares problem
        // Can just directly call a weighted least squares solver, but there's enough
        // structure in our problem it's a little faster just to build the normal eqn
        // matrices ourselves

        double A00 = 0, A01 = 0, A10 = 0, A11 = 0;
        double b0 = 0, b1 = 0;
        double x0 = 0, x1 = 0;

        for (int32_t ei = 0; ei < num_alignments; ++ei) {
            event_alignment_t ea = alignment_output[ei];
            if (ea.hmm_state == 'M') {
                //std::string model_kmer = ea.rc ? pore_model.pmalphabet->reverse_complement(ea.ref_kmer) : ea.ref_kmer;
                uint32_t rank = get_kmer_rank(ea.ref_kmer, KMER_SIZE);

                double raw_event = et.event[ea.event_idx].mean;
                double level_mean = pore_model[rank].level_mean;
                double level_stdv = pore_model[rank].level_stdv;

                double inv_var = 1. / (level_stdv * level_stdv);
                double mu = level_mean;
                double e = raw_event;

                A00 += inv_var;
                A01 += mu * inv_var;
                A11 += mu * mu * inv_var;

                b0 += e * inv_var;
                b1 += mu * e * inv_var;
            }
        }

        A10 = A01;

        // perform the linear solve
        //Eigen::VectorXd x = A.fullPivLu().solve(b);
        double div = A00 * A11 - A01 * A10;
        x0 = -(A01 * b1 - A11 * b0) / div;
        x1 = (A00 * b1 - A10 * b0) / div;

        double shift = x0;
        double scale = x1;
        //double drift = 0.;
        double var = 1.0;

        if (scale_var) {
            var = 0.;
            for (int32_t ei = 0; ei < num_alignments; ++ei) {
                event_alignment_t ea = alignment_output[ei];
                if (ea.hmm_state == 'M') {
                    uint32_t rank = get_kmer_rank(ea.ref_kmer, KMER_SIZE);
                    double raw_event = et.event[ea.event_idx].mean;
                    double level_mean = pore_model[rank].level_mean;
                    double level_stdv = pore_model[rank].level_stdv;
                    double yi = (raw_event - shift - scale * level_mean);
                    var += yi * yi / (level_stdv * level_stdv);
                }
            }
            var /= num_M_state;
            var = sqrt(var);
        }

        scallings->shift = shift;
        scallings->scale = scale;
        //scallings->drift=drift;
        scallings->var = var;
#ifdef CACHED_LOG
        scallings->log_var = log(var);
#endif

        recalibrated = true;

#ifdef DEBUG_RECALIB_SCALING
        fprintf(stderr, "shift: %.2lf scale: %.2lf var: %.2lf\n",
                scallings->shift, scallings->scale, scallings->var);
        //fprintf(stderr, "truth shift: %.2lf scale: %.2lf\n", pore_model.shift, pore_model.scale);
#endif
    }

    return recalibrated;
}
