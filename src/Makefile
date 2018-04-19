#OPTIMIZE = -pg -Wall
OPTIMIZE = -O3 -g0
#OPTIMIZE = -O2 -g0
CC = gcc

elk : coalition.o characteristic.o dp.o idp.o grasp.o ip.o partitions.o elk.o sandholm.o time_utils.o
	${CC} ${OPTIMIZE} coalition.o characteristic.o dp.o idp.o ip.o partitions.o grasp.o sandholm.o time_utils.o elk.o -lm -o elk
elk.o: elk.c
	${CC} ${OPTIMIZE} -c elk.c
coalition.o : coalition.c coalition.h
	${CC} ${OPTIMIZE} -c coalition.c
characteristic.o : characteristic.c  characteristic.h
	${CC} ${OPTIMIZE} -c characteristic.c
grasp.o : grasp.c grasp.h
	${CC} ${OPTIMIZE} -c grasp.c
dp.o : dp.c dp.h
	${CC} ${OPTIMIZE} -c dp.c
idp.o : idp.c idp.h
	${CC} ${OPTIMIZE} -c idp.c
ip.o : ip.c ip.h
	${CC} ${OPTIMIZE} -c ip.c
partitions.o : partitions.c partitions.h
	${CC} ${OPTIMIZE} -c partitions.c
sandholm.o : sandholm.c sandholm.h
	${CC} ${OPTIMIZE} -c sandholm.c
time_utils.o : time_utils.c time_utils.h
	${CC} ${OPTIMIZE} -c time_utils.c

clean:
	rm -f elk *.o 
