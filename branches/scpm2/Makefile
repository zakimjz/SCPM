CC = g++

CFLAGS = -Wall -g -lm -lpthread -O3 #-pg

RM = rm -f

SCPM_OBJS = main.o str_corr.o quasi_clique.o subgraph.o  dictionary.o graph.o io.o util.o getopt_pp.o
NAIVE_OBJS = naive.o str_corr.o quasi_clique.o subgraph.o  dictionary.o graph.o io.o util.o getopt_pp.o
EXP_EPS_OBJS = sim_analyt_str_corr.o str_corr.o quasi_clique.o subgraph.o  dictionary.o graph.o io.o util.o getopt_pp.o

all: scpm naive sim_analyt_str_corr

scpm: $(SCPM_OBJS)
	@echo ""
	@echo " --- scpm ---"
	@$(CC) $(CFLAGS) $(SCPM_OBJS) -o scpm 
	@echo ""

naive: $(NAIVE_OBJS)
	@echo ""
	@echo " --- naive ---"
	@$(CC) $(CFLAGS) $(NAIVE_OBJS) -o naive
	@echo ""

sim_analyt_str_corr: $(EXP_EPS_OBJS)
	@echo ""
	@echo " --- sim_analyt_str_corr ---"
	@$(CC) $(CFLAGS) $(EXP_EPS_OBJS) -o sim_analyt_str_corr
	@echo ""

%.o: %.cc %.h
	@echo " --- COMPILING OBJECT \"$@\""
	@$(CC) $(CFLAGS) $< -c  

clean:
	$(RM) scpm naive sim_analyt_str_corr *.o
	clear

