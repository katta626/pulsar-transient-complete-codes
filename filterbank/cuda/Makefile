NVCC ?= nvcc
ARCH ?= -arch=native
FLAGS = -O3 -std=c++17 -Xcompiler -fPIC -Iinclude $(ARCH)

build/libfbcuda.so: src/fbcuda.cu include/fbcuda.h
	mkdir -p build
	$(NVCC) $(FLAGS) -shared -o $@ src/fbcuda.cu -lcufft

clean:
	rm -rf build
