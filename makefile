CC = mcc
FLAGS = -R -nodisplay -R -singleCompThread -R -nojvm

SRC = code
S = $(SRC)/simulate
A = $(SRC)/analysis
B = bin

main: $B/Bridge_Pre_Allocation $B/SimulateDay
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

lightClean:
	rm -f mccExcludedFiles.log  requiredMCRProducts.txt readme.txt

clean:
	rm -rf bin/*
