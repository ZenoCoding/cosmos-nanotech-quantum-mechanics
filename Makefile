CC ?= gcc
PYTHON ?= python3

CPPFLAGS := -D_DEFAULT_SOURCE -D_POSIX_C_SOURCE=200809L
CFLAGS := -std=c11 -O2 -Wall -Wextra
MATH_LIBS := -lm
QUANTUM_LIBS := -llapack -lblas -lm

LAB12_SRC := labs/12/src
LAB12_BUILD := build/12

EARLY_LABS := 00 01 02 03 04 05 06 07 08 09 10 11
LAB00_PROGRAMS := hello snake
LAB01_PROGRAMS := add count decision1 mycount
LAB02_PROGRAMS := adder2 ex newton1 quadratic1
LAB03_PROGRAMS := gravity gravity1
LAB04_PROGRAMS := mass_spring nonideal_spring
LAB05_PROGRAMS := projectilemotion
LAB06_PROGRAMS := kepler kepler2
LAB07_PROGRAMS := 3body
LAB08_PROGRAMS := diffusion diffusion2
LAB09_PROGRAMS := random randomwalk
LAB10_PROGRAMS := dotproduct matrix matrixtimed
LAB11_PROGRAMS := solid

LAB12_UTILITY_PROGRAMS := \
	annealbarrier annealhex evolveannealed evolvebarrier \
	timeevolve timeevolve2d timehex

LAB12_MATRIX_PROGRAMS := tight tight2 tight2d

EARLY_TARGETS :=

define EARLY_LAB_RULES
LAB$(1)_TARGETS := $$(addprefix build/$(1)/,$$(LAB$(1)_PROGRAMS))
EARLY_TARGETS += $$(LAB$(1)_TARGETS)

$$(LAB$(1)_TARGETS): build/$(1)/%: labs/$(1)/src/%.c | build/$(1)
	$$(CC) $$(CPPFLAGS) $$(CFLAGS) $$< $$(MATH_LIBS) -o $$@

.PHONY: lab$(1)
lab$(1): $$(LAB$(1)_TARGETS)
endef

$(foreach lab,$(EARLY_LABS),$(eval $(call EARLY_LAB_RULES,$(lab))))

LAB12_UTILITY_TARGETS := $(addprefix $(LAB12_BUILD)/,$(LAB12_UTILITY_PROGRAMS))
LAB12_MATRIX_TARGETS := $(addprefix $(LAB12_BUILD)/,$(LAB12_MATRIX_PROGRAMS))
ALL_TARGETS := $(EARLY_TARGETS) $(LAB12_UTILITY_TARGETS) $(LAB12_MATRIX_TARGETS)

.PHONY: all early lab12 visuals check clean

all: $(ALL_TARGETS)

early: $(EARLY_TARGETS)

lab12: $(LAB12_UTILITY_TARGETS) $(LAB12_MATRIX_TARGETS)

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

build/%:
	mkdir -p $@

visuals:
	MPLCONFIGDIR=/tmp/cosmos-matplotlib $(PYTHON) scripts/generate_gallery.py

check: all
	$(PYTHON) -m compileall -q labs scripts
	@printf 'Built %s C programs successfully.\n' "$(words $(ALL_TARGETS))"

clean:
	rm -rf build
