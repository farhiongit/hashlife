TEMPLATES_DIR = ../templates
#CC=clang -fmacro-backtrace-limit=0 -Wno-gnu-binary-literal
CC=gcc
CFLAGS += -Wall -Wextra -Werror -I$(TEMPLATES_DIR) -Wno-format -Wno-format-security
#CFLAGS += -pedantic
#CFLAGS += -fdebug-cpp
#LDFLAGS += -pthread
CFLAGS += -O3
#CFLAGS += -g -DDEBUG -DTU
#CFLAGS += -pg
#LDFLAGS += -pg

all: bitl.o bitl hgolbi.o hgolbi_example infos

.PHONY: infos
infos:
	cloc --by-file hgolbi.h hgolbi.c bitl.h bitl.c
	nm -g --defined-only ./hgolbi.o
	nm -g --defined-only ./bitl.o
	time ./hgolbi_example -U -x-9_10,3_4 -y-5_6,7_8 -t1_0 -t2_0 </dev/null

bitl.o: bitl.c bitl.h

bitl: bitl.c
	gcc -Wno-format -Wall -Wextra -Werror -I$(TEMPLATES_DIR) -g -DDEBUG -DTUEXE -o bitl bitl.c

#hgolbi.o: CFLAGS += -DDEBUG
hgolbi.o: hgolbi.h hgolbi.c bitl.o $(TEMPLATES_DIR)/set_impl.h $(TEMPLATES_DIR)/bnode_impl.h $(TEMPLATES_DIR)/bnode.h $(TEMPLATES_DIR)/vfunc.h $(TEMPLATES_DIR)/set.h $(TEMPLATES_DIR)/defops.h $(TEMPLATES_DIR)/list_impl.h $(TEMPLATES_DIR)/list.h

hgolbi_example.o: hgolbi_example.c hgolbi.h bitl.h

hgolbi_example: hgolbi.o bitl.o hgolbi_example.o
