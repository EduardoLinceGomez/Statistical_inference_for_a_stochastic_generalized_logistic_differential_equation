# Disable the default rules
MAKEFLAGS += --no-builtin-rules --no-builtin-variables

# Project name
NAME := sde_gen_log
DEBUGG_COMPILE_OPTIONS := -g -c -i8 -I"${MKLROOT}/include"
COMPILE_OPTIONS := -c -i8 -I"${MKLROOT}/include"
LINK_OPTIONS = -L${MKLROOT}/lib/intel64
LINK_OPTIONS += -lmkl_intel_ilp64 -lmkl_sequential
LINK_OPTIONS += -lmkl_core -lpthread -lm -ldl

# Configuration settings
# FC := gfortran
# Debugg symbols
FC := ifort
AR := ar rcs
LD := $(FC)
RM := rm -f
GD := ./gen_dep.awk

# List of all source files
SRCS := src/mod_random_number_generator.f90 \
	src/mod_stochastic_logistic_model.f90 \
	src/mod_milstein_solver.f90 \
	src/mod_brownian_motion.f90 \
	src/mod_path_sampler.f90 \
	src/mod_montecarlo_path_sampler.f90

TEST_SRCS := tests/mod_random_number_generator_test.f90 \
	tests/mod_stochastic_logistic_model_test.f90 \
	tests/mod_milstein_solver_test.f90 \
	tests/mod_milstein_sampler_test.f90 \
	tests/mod_path_sampler_test.f90 \
	tests/mod_brownian_motion_test.f90 \
	tests/mod_montecarlo_path_sampler_test.f90
# Add source and tests directories to search paths
vpath % .: src
vpath % .: tests

# Define a map from each file name to its object file
obj = $(src).o
$(foreach src, $(SRCS) $(TEST_SRCS), $(eval $(src) := $(obj)))

# Create lists of the build artefacts in this project
OBJS := $(addsuffix .o, $(SRCS))
DEPS := $(addsuffix .d, $(SRCS))
TEST_OBJS := $(addsuffix .o, $(TEST_SRCS))
TEST_DEPS := $(addsuffix .d, $(TEST_SRCS))
LIB := $(patsubst %, lib%.a, $(NAME))
TEST_BIN := $(patsubst %.f90, %.out, $(TEST_SRCS))

.PHONY: all clean
all: $(LIB) $(TEST_BIN)
# all: $(TEST_BIN)

# Create the static library from the object files
$(LIB): $(OBJS)
	$(AR) $@ $^

$(TEST_BIN): %.out: %.f90.o $(LIB)
	$(LD) $(LINK_OPTIONS) -o $@ $^

# Create object files from Fortran source
$(OBJS) $(TEST_OBJS): %.o: % | %.d
	$(FC) $(DEBUGG_COMPILE_OPTIONS) -o $@ $<

# Process the Fortran source for module dependencies
$(DEPS) $(TEST_DEPS): %.d: %
	$(GD) $< > $@

# Define all module interdependencies
include $(DEPS) $(TEST_DEPS)
$(foreach dep, $(OBJS) $(TEST_OBJS), $(eval $(dep): $($(dep))))

# Cleanup, filter to avoid removing source code by accident
make run:
	make
	./tests/mod_random_number_generator_test.out \
	./tests/mod_stochastic_logistic_model_test.out \
	./tests/mod_milstein_solver_test.out \
	./tests/mod_milstein_sampler_test.out \
	./tests/mod_path_sampler_test.out

clean:
	$(RM) $(filter %.o, $(OBJS) $(TEST_OBJS)) \
	$(filter %.d, $(DEPS) $(TEST_DEPS)) \
	$(filter %.out, $(TEST_BIN)) $(LIB)#\
	$(wildcard *.mod)