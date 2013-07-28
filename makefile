all:
	#gcc -O3 -fopenmp natestm.cpp graphics.cpp -o test -lGL -lGLU -lglut -lglui -lm -lstdc++
	#gcc -O3 natestm.cpp graphics.cpp -o test -lGL -lGLU -lglut -lglui -lm -lstdc++
	nvcc -O3 natestm.cpp graphics.cpp g_kernel.cu -o test -lGL -lGLU -lglut -lglui -lm -lstdc++