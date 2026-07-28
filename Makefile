CC ?= gcc
PYTHON ?= python3

CPPFLAGS := -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L
CFLAGS := -std=c11 -O2 -Wall -Wextra
MATH_LIBS := -lm
QUANTUM_LIBS := -llapack -lblas -lm

EARLY_SRC := labs/01-11/src
LAB12_SRC := labs/12/src
EARLY_BUILD := build/01-11
LAB12_BUILD := build/12

EARLY_PROGRAMS := \
	3body add adder2 count decision1 diffusion diffusion2 dotproduct ex \
	gravity gravity1 hello kepler kepler2 mass_spring matrix matrixtimed \
	mycount newton1 nonideal_spring projectilemotion quadratic1 random \
	randomwalk snake solid

LAB12_UTILITY_PROGRAMS := \
	annealbarrier annealhex evolveannealed evolvebarrier \
	timeevolve timeevolve2d timehex

LAB12_MATRIX_PROGRAMS := tight tight2 tight2d

EARLY_TARGETS := $(addprefix $(EARLY_BUILD)/,$(EARLY_PROGRAMS))
LAB12_UTILITY_TARGETS := $(addprefix $(LAB12_BUILD)/,$(LAB12_UTILITY_PROGRAMS))
LAB12_MATRIX_TARGETS := $(addprefix $(LAB12_BUILD)/,$(LAB12_MATRIX_PROGRAMS))
ALL_TARGETS := $(EARLY_TARGETS) $(LAB12_UTILITY_TARGETS) $(LAB12_MATRIX_TARGETS)

.PHONY: all early lab12 visuals check clean

all: $(ALL_TARGETS)

early: $(EARLY_TARGETS)

lab12: $(LAB12_UTILITY_TARGETS) $(LAB12_MATRIX_TARGETS)

$(EARLY_TARGETS): $(EARLY_BUILD)/%: $(EARLY_SRC)/%.c | $(EARLY_BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(MATH_LIBS) -o $@

$(LAB12_UTILITY_TARGETS): $(LAB12_BUILD)/%: \
		$(LAB12_SRC)/%.c $(LAB12_SRC)/utilities.c $(LAB12_SRC)/utilities.h \
		| $(LAB12_BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LAB12_SRC)/utilities.c \
		$(QUANTUM_LIBS) -o $@

$(LAB12_MATRIX_TARGETS): $(LAB12_BUILD)/%: \
		$(LAB12_SRC)/%.c $(LAB12_SRC)/matrix.c $(LAB12_SRC)/matrix.h \
		$(LAB12_SRC)/utilities.c $(LAB12_SRC)/utilities.h | $(LAB12_BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LAB12_SRC)/matrix.c \
		$(LAB12_SRC)/utilities.c $(QUANTUM_LIBS) -o $@

$(EARLY_BUILD) $(LAB12_BUILD):
	mkdir -p $@

visuals:
	MPLCONFIGDIR=/tmp/cosmos-matplotlib $(PYTHON) scripts/generate_gallery.py

check: all
	$(PYTHON) -m compileall -q labs scripts
	@printf 'Built %s C programs successfully.\n' "$(words $(ALL_TARGETS))"

clean:
	rm -rf build
