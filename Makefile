all:
	gcc studentsseq.c -o studentsseq -lm -fopenmp
	./studentsseq < tests/1.in
	gcc studentspar.c -o studentspar -lm -fopenmp
	./studentspar < tests/1.in