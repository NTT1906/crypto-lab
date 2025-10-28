#include <stdio.h>
#include <stdint.h>
#include <time.h>

typedef unsigned int uint;
#define SIZEOF_UINT (sizeof(uint) * 8)

#define N 200000000

// Branch version
static inline uint u1_shr_branch(uint num, uint amnt) {
	if (amnt >= SIZEOF_UINT) return 0;
	return num >> amnt;
}

// Branchless version
static inline uint u1_shr_bl(uint num, uint amnt) {
	return (num >> (amnt & (SIZEOF_UINT - 1))) & -(amnt < SIZEOF_UINT);
}

int main() {
	uint sum = 0;
	clock_t start, end;

	// Test data
	uint nums[4] = {0xFFFFFFFF, 0x12345678, 0x0F0F0F0F, 0x80000000};
	uint shifts[4] = {0, 5, 31, 40}; // include out-of-range

	// Branch version
	start = clock();
	for (int i = 0; i < N; i++) {
		sum += u1_shr_branch(nums[i & 3], shifts[i & 3]);
	}
	end = clock();
	printf("Branch version: %f sec, sum=%u\n", (double)(end - start)/CLOCKS_PER_SEC, sum);

	sum = 0;
	// Branchless version
	start = clock();
	for (int i = 0; i < N; i++) {
		sum += u1_shr_bl(nums[i & 3], shifts[i & 3]);
	}
	end = clock();
	printf("Branchless version: %f sec, sum=%u\n", (double)(end - start)/CLOCKS_PER_SEC, sum);

	return 0;
}
