all:
	gcc studentsseq.c -o studentsseq -lm
	./studentsseq < tests/1.in
	gcc studentspar.c -o studentspar -lm -fopenmp
	./studentspar < tests/1.in