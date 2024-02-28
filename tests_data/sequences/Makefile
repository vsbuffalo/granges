# does all the bgzip/indexing of test files

.PHONY: clean all

all: test_case_01.fa.gz.fai test_case_01.fa.gz.gzi


test_case_01.fa.gz: test_case_01.fa
	bgzip -c $< > $@


test_case_01.fa.gz.fai test_case_01.fa.gz.gzi: test_case_01.fa.gz
	samtools faidx $< > $@

clean: 
	rm -f test_case_01.fa.gz.fai test_case_01.fa.gz.gzi test_case_01.fa.gz.fai test_case_01.fa.gz


