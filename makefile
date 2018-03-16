CC = mcc
FLAGS = -R -nodisplay -R -singleCompThread -R -nojvm

SRC = code
SIM = $(SRC)/simulate
ANA = $(SRC)/analysis
BIN = bin

main: $(BIN)/Bridge_Pre_Allocation $(BIN)/SimulateDay $(BIN)/createDataMatrix $(BIN)/zero_one_scale
	rm -f mccExcludedFiles.log
	rm -f requiredMCRProducts.txt
	rm -f readme.txt

$(BIN)/SimulateDay: $(SIM)/SimulateDay.m
	$(CC) $(FLAGS) -m $< -o SimulateDay
	rm -f run_SimulateDay.sh
	mv SimulateDay $@

$(BIN)/Bridge_Pre_Allocation: $(SIM)/Bridge_Pre_Allocation.m
	$(CC) $(FLAGS) -m $< -o Bridge_Pre_Allocation
	rm -f run_Bridge_Pre_Allocation.sh
	mv Bridge_Pre_Allocation $@

$(BIN)/createDataMatrix: $(ANA)/createDataMatrix.m $(ANA)/day2Mat.m
	$(CC) $(FLAGS) -m $< -o createDataMatrix
	rm -f run_createDataMatrix.sh
	mv createDataMatrix $@

$(BIN)/zero_one_scale: $(ANA)/zero_one_scale.m
	$(CC) -R -nodisplay -R -nojvm -m $< -g -o zero_one_scale
	rm -f run_zero_one_scale.sh
	mv zero_one_scale $@

lightClean:
	rm -f mccExcludedFiles.log  requiredMCRProducts.txt readme.txt

clean:
	rm -rf bin/*
