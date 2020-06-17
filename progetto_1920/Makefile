main_iplib: main_iplib.c ip_lib.o bmp.o
	gcc main_iplib.c ip_lib.o bmp.o -omain_iplib -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

ip_lib.o: ip_lib.c ip_lib.h bmp.o
	gcc ip_lib.c -oip_lib.o --ansi --pedantic -Wall -c -lm

bmp.o: bmp.c bmp.h
	gcc bmp.c -obmp.o -Wall -ggdb -c -lm

clean:
	@rm -f main_iplib.o ip_lib.o bmp.o 
