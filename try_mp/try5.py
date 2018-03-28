import 

print(gmpy2.get_context())

gmpy2.get_context(precision=200)

print(gmpy2.sqrt(gmpy2.mpc("1+2j")))