test_bmp: test_bmp.o bmp.o
	gcc test_bmp.o bmp.o -otest_bmp -lm

test_bmp.o: test_bmp.c bmp.h
	gcc test_bmp.c -otest_bmp.o -Wall -c -lm

bmp.o: bmp.c bmp.h
	gcc bmp.c -obmp.o -Wall -c -lm

clean:
	@rm -f bmp.o test_bmp.o test_bmp
