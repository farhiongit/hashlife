TEMPLATES_DIR = ../templates
#CC=clang -fmacro-backtrace-limit=0 -Wno-gnu-binary-literal
#CC=gcc
CFLAGS += -O3

all: bitl.o hgolbi.o infos TU

.PHONY: infos
infos: bitl.o hgolbi.o
	cloc --by-file hgolbi.h hgolbi.c bitl.h bitl.c
	nm -g --defined-only ./hgolbi.o
	nm -g --defined-only ./bitl.o

.PHONY: TU
TU: bitl hgolbi_example
	./bitl
	time -v ./hgolbi_example -U -x-9_10,3_4 -y-5_6,7_8 -t1_0 -t2_0 </dev/null

bitl.o: CFLAGS += -Wno-format -Wno-format-security
bitl.o: bitl.c bitl.h

bitl: CFLAGS += -Wno-format -Wno-format-security -g -DDEBUG -DTUEXE
bitl: bitl.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

#hgolbi.o: CFLAGS += -DDEBUG
hgolbi.o: CFLAGS += -I$(TEMPLATES_DIR)
hgolbi.o: hgolbi.h hgolbi.c bitl.o $(TEMPLATES_DIR)/set_impl.h $(TEMPLATES_DIR)/bnode_impl.h $(TEMPLATES_DIR)/bnode.h $(TEMPLATES_DIR)/vfunc.h $(TEMPLATES_DIR)/set.h $(TEMPLATES_DIR)/defops.h $(TEMPLATES_DIR)/list_impl.h $(TEMPLATES_DIR)/list.h

hgolbi_example.o: CFLAGS += -Wno-format -Wno-format-security  # since register_printf_specifier (non-standard gnu extension) is used.
hgolbi_example.o: hgolbi_example.c hgolbi.h bitl.h

hgolbi_example: hgolbi.o bitl.o hgolbi_example.o
