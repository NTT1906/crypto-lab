# Note:
Final version: /new/bigint.h

Features:
- add/sub/mul/div/mod bigint
- add_mod, mul_mod, pow_mod (naive), pow_mod (montgomery - 10x naive)
- shift left, shift right, shift left mod
- bui_from_dec(): decimal string to bigint BI_BIT bit (512)
- bui_to_dec(): bigint to decimal string
- bul_to_dec(): long bigint to decimal string
- bul_to_hex(): bigint to hex string
- bul_high, bul_low: MSH and MSH of long bigint
Incomplete:
- bui_from_hex()
- karatsuba multiplication
- divmod_knuth(): Knuth Algorithm D
- poly mul
- poly pow
- aks (in new bigint)