CC = mcc
FLAGS = -R -nodisplay -R -singleCompThread -R -nojvm

SRC = code
S = $(SRC)/simulate
A = $(SRC)/analysis
B = bin

main: $B/Bridge_Pre_Allocation $B/SimulateDay $B/createDataMatrix $B/zero_one_scale
	rm -f mccExcludedFiles.log
	rm -f requiredMCRProducts.txt
	rm -f readme.txt

$B/SimulateDay: $S/SimulateDay.m
	$(CC) $(FLAGS) -m $< -o SimulateDay
	rm -f run_SimulateDay.sh
	mv SimulateDay $@

$B/Bridge_Pre_Allocation: $S/Bridge_Pre_Allocation.m
	$(CC) $(FLAGS) -m $< -o Bridge_Pre_Allocation
	rm -f run_Bridge_Pre_Allocation.sh
	mv Bridge_Pre_Allocation $@

$B/createDataMatrix: $A/createDataMatrix.m $A/day2Mat.m
	$(CC) $(FLAGS) -m $< -o createDataMatrix
	rm -f run_createDataMatrix.sh
	mv createDataMatrix $@

$B/zero_one_scale: $A/zero_one_scale.m
	$(CC) -R -nodisplay -R -nojvm -m $< -g -o zero_one_scale
	rm -f run_zero_one_scale.sh
	mv zero_one_scale $@

lightClean:
	rm -f mccExcludedFiles.log  requiredMCRProducts.txt readme.txt

clean:
	rm -rf bin/*
