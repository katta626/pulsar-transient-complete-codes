VENDOR := ../vendor/psrdada
UNAME := $(shell uname)
EXT := $(if $(filter Darwin,$(UNAME)),dylib,so)
LIB := build/libfbdada.$(EXT)
SRCS := src/fbdada.c $(addprefix $(VENDOR)/,ipcbuf.c ipcio.c ipcutil.c dada_hdu.c ascii_header.c multilog.c tmutil.c)
CFLAGS ?= -O2
CFLAGS += -fPIC -Iinclude -I$(VENDOR) -w

all: $(LIB)

$(LIB): $(SRCS)
	mkdir -p build
	$(CC) $(CFLAGS) -shared -o $@ $(SRCS) -lpthread

clean:
	rm -rf build

.PHONY: all clean
