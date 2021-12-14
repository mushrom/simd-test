#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <functional>

#define NUM (1 << 22)

#if defined(__ARM_ARCH_5__) \
	|| defined(__ARM_ARCH_6__) \
	|| defined(__ARM_ARCH_7__)
#define IS_ARM
#endif

#if defined(IS_ARM)
#warning asdf!
#include <arm_neon.h>
#endif

static inline uint64_t get_timestamp( void ){
#if defined(__x86_64__)
	register uint64_t a, d;
	asm volatile ( "rdtsc" : "=a"(a), "=d"(d));
	return (d << 32) | a;
#else
	return clock();
#endif
}

static inline double get_freq(void) {
#if defined(__x86_64__)
	uint64_t time_a = get_timestamp();
	usleep(1000000);
	//sleep(1);
	//for ( volatile int i = 0; i < 10000000; i++ );
	uint64_t time_b = get_timestamp();
	uint64_t delta = time_b - time_a;

	printf( "total real clocks: %lu\n", delta );
	//printf( "timestamp a: %lu, timestamp b: %lu\n", time_a, time_b );

	double ghertz = delta / 1000000000.0;

	printf( "approx. cpu speed: %.4lfghz\n", ghertz );
	return delta;
#else
	return CLOCKS_PER_SEC;
#endif
}

struct quadratic {
	float a;
	float b;
	float c;
};

struct quadraticSoA {
	float a[NUM];
	float b[NUM];
	float c[NUM];
};

struct quadraticAoSoA {
	float a[8];
	float b[8];
	float c[8];
};

typedef float v2sf __attribute__((vector_size (8)));
typedef float v4sf __attribute__((vector_size (16)));
typedef float v8sf __attribute__((vector_size (32)));

// using the same sqrt instruction for all calculations for consistency
static inline float ssesqrtf(float n) {
#if defined(__i386__)
	v4sf temp = {n}; 
	// loading into 4-wide vector, don't think there's a way to do sqrt
	// on only a single value
	return __builtin_ia32_sqrtps(temp)[0];
#elif defined(IS_ARM)
	v2sf temp = {n}; 
	// use reciprocal sqrt approximation, sqrt is on 64 bit arm only it seems
	return 1.f/vrsqrte_f32(temp)[0];
#else
	return sqrtf(n);
#endif
}

void doPlainAoS(quadratic *quads, float *results) {
	for (int i = 0; i < NUM; i++) {
		quadratic& q = quads[i];
		//results[i] = (-q.b + sqrtf(q.b*q.b - 4*q.a*q.c)) / (2*q.a);
		results[i] = (-q.b + ssesqrtf(q.b*q.b - 4*q.a*q.c)) / (2*q.a);
	}
}

//float *doPlainSoA(quadraticSoA& quads) __attribute__((optimize("-fnotree-loop-vectorize -fnotree-vectorize")));

//void doPlainSoA(float *a, float *b, float *c, float *results) __attribute__((optimize("-fno-tree-vectorize")));
void doPlainSoA(quadraticSoA *quads, float *results) {
	float *a = quads->a;
	float *b = quads->b;
	float *c = quads->c;

	for (int n = 0; n < NUM; n++) {
		//results[n] = (-b[n] + sqrtf(b[n]*b[n] - 4*a[n]*c[n])) / (2*a[n]);
		results[n] = (-b[n] + ssesqrtf(b[n]*b[n] - 4*a[n]*c[n])) / (2*a[n]);
	}
}

#if (defined(__i386__) || defined(__x86_64__)) && defined(USESSE)
void doSoA_sse(quadraticSoA *quads, float *results) {
	quadraticSoA *q = quads;

	for (int n = 0; n < NUM; n += 4) {
		v4sf a = __builtin_ia32_loadups(q->a + n);
		v4sf b = __builtin_ia32_loadups(q->b + n);
		v4sf c = __builtin_ia32_loadups(q->c + n);

		v4sf temp = b*b - 4*a*c;
		v4sf sqr  = __builtin_ia32_sqrtps(temp);
		v4sf result = (-b*sqr) / 2*a;

		__builtin_ia32_storeups(results+n, result);
	}
}

void doAoSoA_sse(quadraticAoSoA *quads, float *results) {
	for (int n = 0; n < NUM; n += 8) {
		quadraticAoSoA *q = quads + n/8;

		for (int k = 0; k < 8; k += 4) {
			v4sf a = __builtin_ia32_loadups(q->a + k);
			v4sf b = __builtin_ia32_loadups(q->b + k);
			v4sf c = __builtin_ia32_loadups(q->c + k);

			v4sf temp = b*b - 4*a*c;
			v4sf sqr  = __builtin_ia32_sqrtps(temp);
			v4sf result = (-b*sqr) / 2*a;

			__builtin_ia32_storeups(results+n + k, result);
		}
	}
}
#endif

#if defined(__x86_64__) && defined(USEAVX)
void doSoA_avx(quadraticSoA *quads, float *results) {
	quadraticSoA *q = quads;

	for (int n = 0; n < NUM; n += 8) {
		v8sf a = __builtin_ia32_loadups256(q->a + n);
		v8sf b = __builtin_ia32_loadups256(q->b + n);
		v8sf c = __builtin_ia32_loadups256(q->c + n);

		v8sf temp = b*b - 4*a*c;
		v8sf sqr  = __builtin_ia32_sqrtps256(temp);
		v8sf result = (-b*sqr) / 2*a;

		__builtin_ia32_storeups256(results+n, result);
	}
}

void doAoSoA_avx(quadraticAoSoA *quads, float *results) {
	for (int n = 0; n < NUM; n += 8) {
		quadraticAoSoA *q = quads + n/8;

		v8sf a = __builtin_ia32_loadups256(q->a);
		v8sf b = __builtin_ia32_loadups256(q->b);
		v8sf c = __builtin_ia32_loadups256(q->c);

		v8sf temp = b*b - 4*a*c;
		v8sf sqr  = __builtin_ia32_sqrtps256(temp);
		v8sf result = (-b*sqr) / 2*a;

		__builtin_ia32_storeups256(results+n, result);
	}
}
#endif

#if defined(IS_ARM) && defined(USENEON)
void doSoA_neon4(quadraticSoA *quads, float *results) {
	quadraticSoA *q = quads;

	for (int n = 0; n < NUM; n += 4) {
		v4sf a = vld1q_f32(q->a + n);
		v4sf b = vld1q_f32(q->b + n);
		v4sf c = vld1q_f32(q->c + n);

		v4sf temp = b*b - 4*a*c;
		v4sf sqr  = 1.f/vrsqrteq_f32(temp);
		v4sf result = (-b*sqr) / 2*a;

		vst1q_f32(results + n, result);
	}
}

void doAoSoA_neon4(quadraticAoSoA *quads, float *results) {
	for (int n = 0; n < NUM; n += 8) {
		quadraticAoSoA *q = quads + n/8;

		for (int k = 0; k < 8; k += 4) {
			v4sf a = vld1q_f32(q->a + k);
			v4sf b = vld1q_f32(q->b + k);
			v4sf c = vld1q_f32(q->c + k);

			v4sf temp = b*b - 4*a*c;
			v4sf sqr  = 1.f/vrsqrteq_f32(temp);
			v4sf result = (-b*sqr) / 2*a;

			vst1q_f32(results+n + k, result);
			//__builtin_ia32_storeups(results+n + k, result);
		}
	}
}

void doSoA_neon2(quadraticSoA *quads, float *results) {
	quadraticSoA *q = quads;

	for (int n = 0; n < NUM; n += 2) {
		v2sf a = vld1_f32(q->a + n);
		v2sf b = vld1_f32(q->b + n);
		v2sf c = vld1_f32(q->c + n);

		v2sf temp = b*b - 4*a*c;
		v2sf sqr  = 1.f/vrsqrte_f32(temp);
		v2sf result = (-b*sqr) / 2*a;

		vst1_f32(results + n, result);
	}
}

void doAoSoA_neon2(quadraticAoSoA *quads, float *results) {
	for (int n = 0; n < NUM; n += 8) {
		quadraticAoSoA *q = quads + n/8;

		for (int k = 0; k < 8; k += 2) {
			v2sf a = vld1_f32(q->a + k);
			v2sf b = vld1_f32(q->b + k);
			v2sf c = vld1_f32(q->c + k);

			v2sf temp = b*b - 4*a*c;
			v2sf sqr  = 1.f/vrsqrte_f32(temp);
			v2sf result = (-b*sqr) / 2*a;

			vst1_f32(results+n + k, result);
			//__builtin_ia32_storeups(results+n + k, result);
		}
	}
}
#endif


float randFloat(void) {
	//return 10*((rand()/RAND_MAX)*2.f - 1.f);
	return 10*((rand()/(float)RAND_MAX)*2.f - 1.f);
	//return 10*((prng()/(float)UINT_MAX)*2.f - 1.f);
}

int main(int argc, char *argv[]) {
	quadratic *inputs = new quadratic[NUM];
	quadraticSoA *soainputs = new quadraticSoA;
	quadraticAoSoA *aosoainputs = new quadraticAoSoA[NUM/8];
	float *results = new float[NUM];
	//static quadratic inputs[NUM];
	//static quadraticSoA soainputs;

	printf("finding roots of %u quadratics (%lu MiB of data)\n", NUM, (sizeof(float[3]) * NUM)/1024/1024);

	srand(time(NULL));
	double freq = get_freq();
	double msfreq = freq / 1000;
	printf("Time resolution: 1/%g\n", freq);
	//double freq = get_freq();

	printf("Generating data...\n");
	for (int i = 0; i < NUM; i++) {
		quadratic& q = inputs[i];

		q.a = randFloat();
		q.b = randFloat();
		q.c = randFloat();

		soainputs->a[i] = randFloat();
		soainputs->b[i] = randFloat();
		soainputs->c[i] = randFloat();

		results[i] = 0;
	}

	for (int i = 0; i < NUM/8; i ++) {
		quadraticAoSoA& q = aosoainputs[i];

		for (int k = 0; k < 8; k++) {
			q.a[k] = randFloat();
			q.b[k] = randFloat();
			q.c[k] = randFloat();
		}
	}

	printf("Done, starting\n");

	auto doBenchmark = [&](const char *name, std::function<void()> proc) {
		uint64_t start = get_timestamp();
		proc();
		uint64_t end = get_timestamp();
		uint64_t diff = end - start;
		printf("%24s : %12llu cycles (%g milliseconds)\n", name, diff, diff / msfreq);
	};

	doBenchmark("AoS (no simd)", [&]() { doPlainAoS(inputs, results); });
	doBenchmark("SoA (no simd)", [&]() { doPlainSoA(soainputs, results); });

#if (defined(__i386__) || defined(__x86_64__)) && defined(USESSE)
	doBenchmark(  "SoA (sse)", [&]() {   doSoA_sse(soainputs, results); });
	doBenchmark("AoSoA (sse)", [&]() { doAoSoA_sse(aosoainputs, results); });
#endif

#if defined(__x86_64__) && defined(USEAVX)
	doBenchmark(  "SoA (avx)", [&]() {   doSoA_avx(soainputs, results); });
	doBenchmark("AoSoA (avx)", [&]() { doAoSoA_avx(aosoainputs, results); });
#endif

#if defined(IS_ARM) && defined(USENEON)
	doBenchmark(  "SoA (neon x2)", [&]() {   doSoA_neon2(soainputs, results); });
	doBenchmark("AoSoA (neon x2)", [&]() { doAoSoA_neon2(aosoainputs, results); });
	doBenchmark(  "SoA (neon x4)", [&]() {   doSoA_neon4(soainputs, results); });
	doBenchmark("AoSoA (neon x4)", [&]() { doAoSoA_neon4(aosoainputs, results); });
#endif

	return 0;
}
