#
# Montgomery reduction algorithm (Python)
#
# Copyright (c) 2025 Project Nayuki
# All rights reserved. Contact Nayuki for licensing.
# https://www.nayuki.io/page/montgomery-reduction-algorithm
#

import math, random, unittest

class MontgomeryReducerTest(unittest.TestCase):

    def test_basic(self) -> None:
        for _ in range(3000):
            bitlen: int = random.randint(2, 100)
            mod: int = random.randrange(1 << bitlen, 2 << bitlen) | 1  # Force it to be odd
            mr: MontgomeryReducer = MontgomeryReducer(mod)

            for _ in range(100):
                x: int = random.randrange(mod)
                y: int = random.randrange(mod)
                u: int = mr.convert_in(x)
                v: int = mr.convert_in(y)
                w: int = mr.multiply(u, v)
                if mr.convert_out(w) != x * y % mod:
                    raise AssertionError()

            for _ in range(10):
                x = random.randrange(mod)
                y = random.randrange(mod)
                u = mr.convert_in(x)
                v = mr.pow(u, y)
                if mr.convert_out(v) != pow(x, y, mod):
                    raise AssertionError()



class MontgomeryReducer:
    modulus: int
    reducerbits: int
    reducer: int
    mask: int
    reciprocal: int
    factor: int
    convertedone: int

    def __init__(self, mod: int):
        # Modulus
        if mod < 3 or mod % 2 == 0:
            raise ValueError("Modulus must be an odd number at least 3")
        self.modulus = mod

        # Reducer
        self.reducerbits = (mod.bit_length() // 32 + 1) * 32  # This is a multiple of 8
        if self.reducerbits > 512:
            self.reducerbits = 512
        self.reducer = 1 << self.reducerbits  # This is a power of 256
        self.mask = self.reducer - 1
        assert (self.reducer > mod) and (math.gcd(self.reducer, mod) == 1)

        # Other computed numbers
        self.reciprocal = pow(self.reducer, -1, mod)
        print(self.reciprocal)
        self.factor = (self.reducer * self.reciprocal - 1) // mod
        factor_fake = self.mask // mod
        print(f"Fake factor: {factor_fake}")
        self.convertedone = self.reducer % mod


    # The range of x is unlimited
    def convert_in(self, x: int) -> int:
        # print(f"SL: {x << self.reducerbits}")
        # print(f"TM: {(x << self.reducerbits) % self.modulus}")
        return (x << self.reducerbits) % self.modulus


    # The range of x is unlimited
    def convert_out(self, x: int) -> int:
        return (x * self.reciprocal) % self.modulus


    # Inputs and output are in Montgomery form and in the range [0, modulus)
    def multiply(self, x: int, y: int) -> int:
        mod: int = self.modulus
        assert (0 <= x < mod) and (0 <= y < mod)
        product: int = x * y
        # print(f"x      = {x}")
        # print(f"y      = {y}")
        # print(f"1: p   = {product}")

        temp: int = ((product & self.mask) * self.factor) & self.mask
        # t1: int = product & self.mask
        # print(f"2.1: t1= {t1}")
        # t2: int = t1 * self.factor
        # print(f"2.2: t2= {t2}")
        # temp: int = t2 & self.mask
        # print(f"2: temp= {temp}")

        reduced: int = (product + temp * mod) >> self.reducerbits
        # print(f"3.1: r1= {temp * mod}")
        # print(f"3.2: r2= {product + temp * mod}")
        # print(f"3: redu= {reduced}")
        result: int = reduced if (reduced < mod) else (reduced - mod)
        # print(f"4: resu= {result}")
        assert 0 <= result < mod
        return result


    # Input x (base) and output (power) are in Montgomery form and in the range [0, modulus); input y (exponent) is in standard form
    def pow(self, x: int, y: int) -> int:
        assert 0 <= x < self.modulus
        if y < 0:
            raise ValueError("Negative exponent")
        z: int = self.convertedone
        # print(f"Highest bit of e: {y.bit_length()}\n")
        while y != 0:
            if y & 1 != 0:
                z = self.multiply(z, x)
            x = self.multiply(x, x)
            y >>= 1
        return z

if __name__ == "__main__":
    # (123456789012345678901234567890 ** 65537) % 1000000000000000000000000000037
# 73399183678185556697836825553
    u = 123456789012345678901234567890
    v = 65537
    m = 1000000000000000000000000000037
    # u = 260428835329122752520818469321216072583938198616075453742527759001901820374664228839496959095353854544158481165265966459
    # v = 1312312317639123213
    mr = MontgomeryReducer(m)
    # mr = MontgomeryReducer(8274904334290417405341624571932224150456224549917673444239237760272785701939526927698156030175801211624849856326839256526253153336777911614668501375751381)
    # mr = MontgomeryReducer(6274904334290417405341624571932224150456224549917673444239237760272785701939526927698156030175801211624849856326839256526253153336777911614668501375751381)
    # u = 260428835329122752520818469321216072583938198616075453742527
    # v = 1312312317639123213
    # mr = MontgomeryReducer(108579795932485217312615519053)
    print(f"reducerBits= {mr.reducerbits}")
    print(f"reciprocal = {mr.reciprocal}")
    print(f"mask       = {mr.mask}")
    print(f"reducer    = {mr.reducer}")
    print(f"factor     = {mr.factor}")
    print(f"c1         = {mr.convertedone}")

    cu = mr.convert_in(u)
    # print(f"cu={cu}")
    r = mr.pow(cu, v)
    cr = mr.convert_out(r)
    print(f"r ={r}\ncr={cr}")
    # unittest.main()