#
# Compiler flags
#
CC     = gcc
CFLAGS = -Wall -Wextra -pthread

#
# Project files
#

SRCS = main.c rsa.c bignum.c
OBJS = $(SRCS:.c=.o)
EXE  = main
BENCHSRCS = benchmark.c rsa.c bignum.c
BENCHOBJS = $(BENCHSRCS:.c=.o)
BENCHEXE  = benchmark

#
# Debug build settings
#
DBGDIR = debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -g -O0 -DDEBUG

#
# Release build settings
#
RELDIR = release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELBENCHEXE = $(RELDIR)/$(BENCHEXE)
RELBENCHOBJS = $(addprefix $(RELDIR)/, $(BENCHOBJS))
RELCFLAGS = -O3 -DNDEBUG

.PHONY: all benchmark clean debug prep release remake

# Default build
all: release

#
# Debug rules
#
debug: prep $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CC) $(CFLAGS) $(DBGCFLAGS) -o $(DBGEXE) $^

$(DBGDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(DBGCFLAGS) -o $@ $<

#
# Release rules
#
release: prep $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(RELEXE) $^

benchmark: prep $(RELBENCHEXE)

$(RELBENCHEXE): $(RELBENCHOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) -o $(RELBENCHEXE) $^

$(RELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(RELCFLAGS) -o $@ $<

#
# Other rules
#
prep:
	@mkdir -p $(DBGDIR) $(RELDIR)

remake: clean all

clean:
	rm -f $(RELEXE) $(RELOBJS) $(RELBENCHEXE) $(RELBENCHOBJS) $(DBGEXE) $(DBGOBJS)
