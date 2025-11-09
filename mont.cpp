/**
 * ======================================================================
 * CODE ĐẦY ĐỦ: (A^B) % m
 * Bao gồm cả "Đường chậm" (Naive) và "Đường nhanh" (Montgomery)
 * Được xây dựng cho thư viện 'cigint.h' (Big-Endian / MSW-first).
 * ======================================================================
 *
 * HƯỚNG DẪN BIÊN DỊCH:
 * 1. Đặt file này cùng thư mục với `cigint.h`.
 * 2. Đảm bảo CIGINT_N (ví dụ 16) và CIGINT_IMPLEMENTATION được định nghĩa.
 *
 * Ví dụ (g++):
 * g++ -std=c++17 -DCIGINT_N=16 -DCIGINT_IMPLEMENTATION full_code.cpp -o full_code
 *
 * (Nếu N=16, 512-bit, quá trình setup R^2 % m có thể mất vài giây)
 */

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iomanip> // Dùng std::hex, std::setw, std::setfill
#include <cassert> // Dùng assert
#include <random>    // Dùng để tạo số ngẫu nhiên

#define CIGINT_IMPLEMENTATION
#ifndef CIGINT_N
#define CIGINT_N 2
#endif

#include "cigint.h" // File .h của bạn
// -------------------------


//==============================================================================
// PHẦN 1: STRUCT 1024-BIT (2N) VÀ CÁC HÀM PHỤ TRỢ MSW
//==============================================================================

// Struct 1024-bit (2N) - MSW-first (Big-Endian)
// data[0] = MSW
// data[2N - 1] = LSW
typedef struct Cigint2N {
	u32 data[CIGINT_N * 2];
} Cigint2N;

// Lấy bit từ Cigint2N (MSW-first)
// u32 get_bit_2N_msw(const Cigint2N *a, u32 pos) {
// 	assert(pos < CIGINT_N * 2 * SIZEOF_U32);
// 	size_t data_index = (CIGINT_N * 2) - 1 - (pos / SIZEOF_U32);
// 	return (a->data[data_index] >> (pos % SIZEOF_U32)) & 1u;
// }

u32 get_bit_2N_msw(const Cigint2N *a, u32 pos) {
	assert(pos < CIGINT_N * 2 * SIZEOF_U32);
	size_t data_index = (CIGINT_N * 2) - 1 - (pos / SIZEOF_U32);
	return (a->data[data_index] >> (pos % SIZEOF_U32)) & 1u;
}

// Set bit cho Cigint (MSW-first) (lấy từ logic cigint.h)
void set_bit_N_msw(Cigint *a, u32 pos, u32 val) {
	if (!a) return;
	assert(pos < CIGINT_N * SIZEOF_U32);
	size_t data_index = CIGINT_N - 1 - (pos / SIZEOF_U32);
	u32 mask = (u32)1 << (pos % SIZEOF_U32);
	a->data[data_index] = (a->data[data_index] & ~mask) | ((val & 1u) ? mask : 0u);
}

// Tìm bit cao nhất của Cigint2N (MSW-first)
u32 highest_order_2N_msw(const Cigint2N *num) {
	for (size_t i = 0; i < CIGINT_N * 2; ++i) { // i=0 là MSW
		if (num->data[i] > 0) {
			u32 word = num->data[i];
			u32 pos_in_word = 0;
			for(int bit = 31; bit >= 0; --bit) {
				if((word >> bit) & 1) {
					pos_in_word = bit;
					break;
				}
			}
			return pos_in_word + ((CIGINT_N * 2) - 1 - i) * SIZEOF_U32;
		}
	}
	// Kiểm tra bit 0 (nếu số là 0)
	if (num->data[CIGINT_N * 2 - 1] == 0) return 0;
	return 0; // là 0
}

// Dịch trái 'rem' 1 bit (MSW-first, N-limbs)
void shl_1_N_msw(Cigint* rem) {
	u32 carry = 0;
	for (size_t i = CIGINT_N; i-- > 0;) {
		u32 v = rem->data[i];
		rem->data[i] = (v << 1) | carry;
		carry = v >> 31;
	}
}


//==============================================================================
// PHẦN 2: "ĐƯỜNG CHẬM" (BOOTSTRAP - NHÂN/CHIA 1024-BIT)
//==============================================================================

/**
 * Phép nhân Đầy đủ (MSW-first): C(2N) = A(N) * B(N)
 */
void cigint_mul_full_msw(Cigint2N* C_res, const Cigint* A, const Cigint* B) {
	for(unsigned int & i : C_res->data) i = 0;
	for (ssize_t i = CIGINT_N - 1; i >= 0; --i) {
		u64 carry = 0;
		for (ssize_t j = CIGINT_N - 1; j >= 0; --j) {
			size_t k = (size_t)i + (size_t)j + 1;
			u64 temp = (u64)A->data[i] * B->data[j] + C_res->data[k] + carry;
			C_res->data[k] = (u32)temp;
			carry = temp >> 32;
		}
		C_res->data[i] += (u32)carry;
	}
}

/**
 * Phép chia Đầy đủ (MSW-first): (Q, R) = A(2N) / B(N)
 */
void cigint_divmod_full_msw(Cigint *q, Cigint *r, const Cigint2N *lhs, const Cigint *rhs) {
	assert(!cigint_is0_ref(rhs));
	if (q) *q = CIGINT_ZERO();
	if (r) *r = CIGINT_ZERO();
	u32 top = highest_order_2N_msw(lhs);
	if (top == 0 && get_bit_2N_msw(lhs, 0) == 0) return; // lhs là 0

	for (int bit = (int)top; bit >= 0; --bit) {
		shl_1_N_msw(r);
		u32 carry_bit = get_bit_2N_msw(lhs, (u32)bit);
		r->data[CIGINT_N - 1] |= carry_bit;

		if (cigint_cmp_ref(r, rhs) >= 0) {
			cigint_sub_ref(r, rhs);
			set_bit_N_msw(q, (u32)bit, 1u);
		}
	}
}

/**
 * (ĐƯỜNG CHẬM 1) Tính C = (A * B) % m
 */
void cigint_mul_mod_slow_msw(Cigint* C_res, const Cigint* A, const Cigint* B, const Cigint* m) {
	Cigint2N C_full;
	cigint_mul_full_msw(&C_full, A, B);
	cigint_divmod_full_msw(NULL, C_res, &C_full, m);
}

/**
 * (ĐƯỜNG CHẬM 2) Tính C = (A ^ B_exp) % m
 */
void cigint_mod_pow_slow(Cigint* C_res, const Cigint* A, const Cigint* B_exp, const Cigint* m) {
	Cigint Result = cigint_from_u32(1);
	Cigint Base = *A;
	Cigint Exp = *B_exp;

	u32 total_bits = cigint_highest_order_ref(&Exp);
	if (total_bits == 0 && cigint_is0_ref(&Exp)) {
		*C_res = Result;
		return;
	}

	for (u32 bit = 0; bit < total_bits; ++bit) {
		if (cigint_get_bit_ref(&Exp, bit)) { // Lặp từ LSB
			cigint_mul_mod_slow_msw(&Result, &Result, &Base, m);
		}
		cigint_mul_mod_slow_msw(&Base, &Base, &Base, m);
	}
	*C_res = Result;
}


//==============================================================================
// PHẦN 3: "ĐƯỜNG NHANH" (MONTGOMERY ENGINE)
//==============================================================================

/**
 * Tính "số thần kỳ" m_inv = -m[0]^-1 (mod 2^32).
 */
u32 tinh_m_inv(u32 m0) {
	if ((m0 & 1) == 0) {
		throw std::runtime_error("m[0] (LSW) phai la so le!");
	}
	u32 inv = 1;
	for (int i=0; i < 5; ++i) inv = inv * (2 - m0 * inv);
	return (u32)(0 - inv);
}

// /** [Helper] So sánh T(N+1, LSW) >= m(N, MSW) */
// static bool lsw_array_gte_msw_cigint(const u32* T, const Cigint* m) {
// 	if (T[CIGINT_N] != 0) return true;
// 	for (size_t i = 0; i < CIGINT_N; ++i) {
// 		u32 t_val = T[(CIGINT_N - 1) - i];
// 		u32 m_val = m->data[i];
// 		if (t_val > m_val) return true;
// 		if (t_val < m_val) return false;
// 	}
// 	return true;
// }
//
// /** [Helper] Trừ: T(N+1, LSW) = T - m(N, MSW) */
// static void lsw_array_sub_msw_cigint(u32* T, const Cigint* m) {
// 	u64 borrow = 0;
// 	for (size_t i = 0; i < CIGINT_N; ++i) {
// 		u64 m_i = m->data[(CIGINT_N - 1) - i];
// 		u64 temp = (u64)T[i] + 0x100000000 - m_i - borrow;
// 		T[i] = (u32)temp;
// 		borrow = (temp >> 32) == 0 ? 1 : 0;
// 	}
// 	T[CIGINT_N] = T[CIGINT_N] - (u32)borrow;
// }

// So sánh: return true iff T (N+1 words, LSW at T[0]) >= m (Cigint, MSW at m->data[0])
static bool lsw_array_gte_msw_cigint(const u32* T, const Cigint* m) {
	// Nếu phần cao nhất (beyond N words) != 0 thì T >= m
	if (T[CIGINT_N] != 0) return true;

	// So sánh từ MSW -> LSW: T_word_at_pos = T[(N-1)-i], m->data[i]
	for (size_t i = 0; i < CIGINT_N; ++i) {
		u32 t_val = T[(CIGINT_N - 1) - i];
		u32 m_val = m->data[i];
		if (t_val > m_val) return true;
		if (t_val < m_val) return false;
	}
	return true; // equal => >=
}

// Trừ: T = T - m (T LSW-ordered, m MSW-ordered). Assumes T >= m.
static void lsw_array_sub_msw_cigint(u32* T, const Cigint* m) {
	unsigned int borrow = 0;
	for (size_t i = 0; i < CIGINT_N; ++i) {
		// m_i in LSW order
		u64 m_i = (u64)m->data[(CIGINT_N - 1) - i];
		u64 t_i = (u64)T[i];

		u64 subtrahend = m_i + borrow;
		if (t_i < subtrahend) {
			// need borrow from higher word
			T[i] = (u32)( (1ULL << 32) + t_i - subtrahend );
			borrow = 1;
		} else {
			T[i] = (u32)(t_i - subtrahend);
			borrow = 0;
		}
	}
	// propagate final borrow into T[N] (it should be 0 or 1)
	if (borrow) {
		// since we assumed T >= m, this should only decrement T[N] by 1
		T[CIGINT_N] = T[CIGINT_N] - 1;
	}
}

// /**
//  * (ĐƯỜNG NHANH 1) Phép nhân Montgomery: res = (A * B * R^-1) % m
//  */

void cigint_montmul_ref(Cigint *res, const Cigint *A, const Cigint* B, const Cigint* m, u32 m_inv) {
	// u32 T[CIGINT_N + 2] = {0};
	u32 T[CIGINT_N * 2] = {0};
	for (size_t i = 0; i < CIGINT_N; ++i) {
		u64 A_i = A->data[(CIGINT_N - 1) - i]; // Lấy A LSW
		u64 carry_mul = 0;
		for (size_t j = 0; j < CIGINT_N; ++j) {
			u64 B_j = B->data[(CIGINT_N - 1) - j]; // Lấy B LSW
			u64 temp = A_i * B_j + T[j] + carry_mul;
			T[j] = (u32)temp;
			carry_mul = temp >> 32;
		}
		u64 temp_carry = (u64)T[CIGINT_N] + carry_mul;
		T[CIGINT_N] = (u32)temp_carry;
		T[CIGINT_N + 1] = T[CIGINT_N + 1] + (u32)(temp_carry >> 32);

		u32 u = T[0] * m_inv;
		u64 carry_mod = 0;
		for (size_t j = 0; j < CIGINT_N; ++j) {
			u64 m_j = m->data[(CIGINT_N - 1) - j]; // Lấy m LSW
			u64 temp = (u64)u * m_j + T[j] + carry_mod;
			T[j] = (u32)temp;
			carry_mod = temp >> 32;
		}
		temp_carry = (u64)T[CIGINT_N] + carry_mod;
		T[CIGINT_N] = (u32)temp_carry;
		T[CIGINT_N + 1] = T[CIGINT_N + 1] + (u32)(temp_carry >> 32);

		for (size_t j = 0; j <= CIGINT_N; ++j) T[j] = T[j + 1];
		T[CIGINT_N + 1] = 0;
	}

	if (lsw_array_gte_msw_cigint(T, m)) {
		lsw_array_sub_msw_cigint(T, m);
	}
	for (size_t i = 0; i < CIGINT_N; ++i) {
		res->data[i] = T[(CIGINT_N - 1) - i]; // Ghi lại MSW
	}
}

//==============================================================================
// PHẦN 4: HÀM LŨY THỪA MODULAR (CẢ HAI CÁCH)
//==============================================================================
// (Hàm cigint_mod_pow_slow đã ở Phần 2)

/**
 * (ĐƯỜNG NHANH 2) Tính C = (A ^ B_exp) % m (Dùng Montgomery)
 */
void cigint_mod_pow_mont(Cigint* C_res, const Cigint* A, const Cigint* B_exp, const Cigint* m) {
	// === 1. SETUP (Dùng "đường chậm" 1 lần) ===
	if ((m->data[CIGINT_N - 1] & 1) == 0) {
		throw std::runtime_error("mod_pow_mont yeu cau m la so le!");
	}
	u32 m_inv = tinh_m_inv(m->data[CIGINT_N - 1]);

	// Tính hằng số R2 = (R^2) % m
	Cigint R2_mod_m;
	Cigint TWO = cigint_from_u32(2);
	Cigint EXP = cigint_from_u32((u32)(2 * CIGINT_N * SIZEOF_U32));

	// (Đây là "chi phí" setup 1 lần)
	cigint_mod_pow_slow(&R2_mod_m, &TWO, &EXP, m);

	// === 2. CHUYỂN VÀO MIỀN MONTGOMERY ===
	Cigint A_mont;
	cigint_montmul_ref(&A_mont, A, &R2_mod_m, m, m_inv);

	Cigint Result_mont;
	Cigint ONE = cigint_from_u32(1);
	cigint_montmul_ref(&Result_mont, &ONE, &R2_mod_m, m, m_inv);

	// === 3. VÒNG LẶP NHANH (Square-and-Multiply) ===
	Cigint Exp = *B_exp;
	u32 total_bits = cigint_highest_order_ref(&Exp);
	if (total_bits == 0 && cigint_is0_ref(&Exp)) {
		*C_res = ONE;
		return;
	}

	for (u32 bit = 0; bit < total_bits; ++bit) {
		if (cigint_get_bit_ref(&Exp, bit)) { // Lặp từ LSB
			cigint_montmul_ref(&Result_mont, &Result_mont, &A_mont, m, m_inv);
		}
		cigint_montmul_ref(&A_mont, &A_mont, &A_mont, m, m_inv);
	}

	// === 4. CHUYỂN RA KHỎI MIỀN MONTGOMERY ===
	cigint_montmul_ref(C_res, &Result_mont, &ONE, m, m_inv);
}


//==============================================================================
// PHẦN 5: HÀM MAIN (SO SÁNH)
//==============================================================================
#include "benchmark.h"
// observable sink to block DCE
static volatile u32 g_sink = 0;
inline void consume(const Cigint& x) {
	u32 acc = 0;
	for (unsigned int i : x.data) acc ^= i; g_sink ^= acc;
}

template<typename F>
void bench_func_ref(const std::string& name, F f, const Cigint& a, const Cigint& b, const Cigint& m) {
	volatile F fp = f;
	Cigint res;
	bench(name, [&]{
		Cigint x = a;
		Cigint y = b;
		Cigint z = m;
		(*fp)(&res, &x, &y, &z);
		consume(x);
	});
}

static void cigint_randomize(Cigint &x) {
	std::random_device rd; std::mt19937 gen(rd());
	std::uniform_int_distribution<u32> dist(0, UINT32_MAX);
	size_t limbs = 1 + gen() % CIGINT_N;
	for (u32 &i : x.data) i = 0;
	for (size_t i = limbs; i-- > 0;) x.data[i] = dist(gen);
}

int main() {
	// (Giả sử N=16, 512-bit)
	if (CIGINT_N < 2) {
		std::cerr << "CIGINT_N phai >= 2" << '\n';
		return 1;
	}

	// // A = 5
	// Cigint A = cigint_from_u32(5);
	// // B_exp = 3
	// Cigint B_exp = cigint_from_u32(3);
	// // m = 13 (phải là số lẻ)
	// Cigint m = cigint_from_u32(13);

	Cigint A, B_exp, m;
	A = cigint_from_dec("16932033513393476012");
	B_exp = cigint_from_dec("5896503122027280589");
	m = cigint_from_dec("17026143986254644813");
	// cigint_randomize(A);
	// cigint_randomize(B_exp);
	// cigint_randomize(m);
	// cigint_set_bit_ref(&m, 0, 1);
	// A = cigint_mod(A, m);
	// B_exp = cigint_mod(B_exp, m);

	// Đảm bảo A, B < m
	// Cigint Q_tmp;
	// cigint_divmod_full_msw(&Q_tmp, &A, &Cigint2N{A.data}, &m);
	// cigint_divmod_full_msw(&Q_tmp, &B_exp, &Cigint2N{B_exp.data}, &m);

	std::cout << "So sanh ModPow voi CIGINT_N = " << CIGINT_N << '\n';
	std::cout << "A = "; cigint_print10(A); std::cout << '\n';
	std::cout << "B = "; cigint_print10(B_exp); std::cout << '\n';
	std::cout << "M = "; cigint_print10(m); std::cout << '\n';
	std::cout << "----------------------------------" << '\n';

	try {
		Cigint C_slow;
		std::cout << "Dang tinh (Duong cham)..." << '\n';
		cigint_mod_pow_slow(&C_slow, &A, &B_exp, &m);
		std::cout << "  = "; cigint_print10(C_slow); std::cout << '\n';

		Cigint C_fast;
		std::cout << "\nDang tinh (Duong nhanh - Montgomery)..." << '\n';
		cigint_mod_pow_mont(&C_fast, &A, &B_exp, &m);
		std::cout << "  = "; cigint_print10(C_fast); std::cout << '\n';

		std::cout << "----------------------------------" << '\n';
		if (cigint_cmp_ref(&C_slow, &C_fast) == 0) {
			std::cout << "THANH CONG: Ca hai ket qua giong nhau!" << '\n';
		} else {
			std::cout << "THAT BAI: Ket qua khac nhau!" << '\n';
		}
		// bench_func_ref("slow pow", cigint_mod_pow_slow, A, B_exp, m);
		// bench_func_ref("fast pow", cigint_mod_pow_mont, A, B_exp, m);
	} catch (const std::exception& e) {
		std::cerr << "Loi: " << e.what() << '\n';
		return 1;
	}

	return 0;
}
