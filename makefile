CC = mcc
FLAGS = -R -nodisplay -R -singleCompThread -R -nojvm

SRC = code
SIM = $(SRC)/simulate
ANA = $(SRC)/analysis
BIN = bin

main: $(BIN)/Bridge_Pre_Allocation $(BIN)/SimulateDay $(BIN)/createDataMatrix
	rm mccExcludedFiles.log
	rm requiredMCRProducts.txt

$(BIN)/SimulateDay: $(SIM)/SimulateDay.m
	$(CC) $(FLAGS) -m $< -o SimulateDay
	rm run_SimulateDay.sh
	mv SimulateDay $@

$(BIN)/Bridge_Pre_Allocation: $(SIM)/Bridge_Pre_Allocation.m
	$(CC) $(FLAGS) -m $< -o Bridge_Pre_Allocation
	rm run_Bridge_Pre_Allocation.sh
	mv Bridge_Pre_Allocation $@

$(BIN)/createDataMatrix: $(ANA)/createDataMatrix.m
	$(CC) $(FLAGS) -m $< -o createDataMatrix
	rm run_createDataMatrix.sh
	mv createDataMatrix $@
