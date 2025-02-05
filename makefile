

main:
	gcc src/gflib/gflib.c src/gflib/polinom.c src/gflib/gen_polinom.c src/main.c -o main

coder:
	gcc src/gflib/gflib.c src/gflib/polinom.c src/gflib/gen_polinom.c src/coder.c -o coder

clean:
	rm main coder