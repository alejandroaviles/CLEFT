CC = g++
Option = -O3 -Wall

Objs = main.o corr_func.o integral.o q_depend_funcs.o k_depend_funcs.o input.o prog_bar.o

Out_file = cleft.exe

main : $(Objs)
	$(CC) $(Option) $(Objs) -o $(Out_file)

main.o : main.cpp
	$(CC) $(Option) -c -o main.o main.cpp

corr_func.o : corr_func.cpp corr_func.h
	$(CC) $(Option) -c -o corr_func.o corr_func.cpp


integral.o : integral.cpp integral.h
	$(CC) $(Option) -c -o integral.o integral.cpp


q_depend_funcs.o : q_depend_funcs.cpp q_depend_funcs.h
	$(CC) $(Option) -c -o q_depend_funcs.o q_depend_funcs.cpp 

k_depend_funcs.o : k_depend_funcs.cpp k_depend_funcs.h
	$(CC) $(Option) -c -o k_depend_funcs.o k_depend_funcs.cpp

input.o : input.cpp input.h
	$(CC) $(Option) -c -o input.o input.cpp

prog_bar.o : prog_bar.cpp prog_bar.h
	$(CC) $(Option) -c -o prog_bar.o prog_bar.cpp

clean:
	rm *.o *~ cleft.exe
