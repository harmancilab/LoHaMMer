all: LoHaMMer

CC = g++
comp_flags = -c -O3 -Wall
gsl_flags = -lgsl -lgslcblas -lz -lpthread
exec_name = bin/LoHaMMer
LIB_DIR = src

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

%.o: %.c
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/main.o \
${LIB_DIR}/lhmmr_gimp_utils.o \
${LIB_DIR}/lhmmr_gimp_state_reduction_utils.o \
${LIB_DIR}/lhmmr_xlog_math.o \
${LIB_DIR}/lhmmr_linear_math.o \
${LIB_DIR}/lhmmr_genome_sequence_tools.o \
${LIB_DIR}/lhmmr_mapped_read_tools.o \
${LIB_DIR}/lhmmr_config.o \
${LIB_DIR}/lhmmr_ansi_cli.o \
${LIB_DIR}/lhmmr_histogram.o \
${LIB_DIR}/lhmmr_seed_manager.o \
${LIB_DIR}/lhmmr_signal_track_tools.o \
${LIB_DIR}/lhmmr_ansi_thread.o \
${LIB_DIR}/lhmmr_annot_region_tools.o \
${LIB_DIR}/lhmmr_ansi_string.o \
${LIB_DIR}/lhmmr_gff_utils.o \
${LIB_DIR}/lhmmr_utils.o \
${LIB_DIR}/lhmmr_nomenclature.o \
${LIB_DIR}/lhmmr_variation_tools.o \
${LIB_DIR}/lhmmr_imputation_utils.o \
${LIB_DIR}/lhmmr_rng.o \
${LIB_DIR}/lhmmr_nucleotide.o

LoHaMMer: ${objs}
	${CC} -O3 ${gsl_flags} -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 

