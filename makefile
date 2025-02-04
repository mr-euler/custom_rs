
all: coder

coder:
	gcc src/coder.c src/binaryToDecimal.c src/count_digits.c src/buildGF.c -o main

gf:
	gcc src/gflib/gflib.c src/gflib/polinom.c src/gflib/gen_polinom.c src/main.c -o main