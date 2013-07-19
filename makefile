all:
	gcc -O3 -fopenmp natestm.cpp graphics.cpp -o test -lGL -lGLU -lglut -lglui -lm -lstdc++
