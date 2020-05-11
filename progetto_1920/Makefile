main: fake_main.c ip_lib.o bmp.o
	gcc fake_main.c ip_lib.o bmp.o -omain --ansi --pedantic -lm

ip_lib.o: ip_lib.c ip_lib.h bmp.o
	gcc ip_lib.c -oip_lib.o --ansi --pedantic -Wall -c -lm

bmp.o: bmp.c bmp.h
	gcc bmp.c -obmp.o -Wall -c -lm

clean:
	@rm -f main_iplib.o ip_lib.o bmp.o 
