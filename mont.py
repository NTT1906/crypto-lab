import gmpy2

# Define the base, exponent, and modulus
base     = 16932033513393476012
exponent =  5896503122027280589
modulus  = 17026143986254644813

# Perform modular exponentiation
result = gmpy2.powmod(base, exponent, modulus)

print(result) # 13076168503990567874
