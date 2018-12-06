TESTSDIR=../tests/
SODSHOCKDIR=sod_shock_tube_plots/
# the output of savefig puts those single quotes on the filenames 
# need to escape spaces in order to avoid resursive make calling itself 
SODFILE1=\'Low\ Order\ in\ Space\ Low\ order\ in\ Time\ Relativistic\ Sod\ Problem\ 1.jpg\'
SODFILE2=\'Low\ Order\ in\ Space\ Low\ order\ in\ Time\ Relativistic\ Sod\ Problem\ 2.jpg\'
SODFILE3=\'Low\ Order\ in\ Space\ Low\ order\ in\ Time\ Relativistic\ Sod\ Problem\ 3.jpg\'
SODFILE4=\'Low\ Order\ in\ Space\ Low\ order\ in\ Time\ NonRelativistic\ Sod\ Problem.jpg\'

# need to use := for simply scanned variable expansion equivalent to cocnatenation
SODPLOT1:=$(SODSHOCKDIR)$(SODFILE1)
SODPLOT2:=$(SODSHOCKDIR)$(SODFILE2)
SODPLOT3:=$(SODSHOCKDIR)$(SODFILE3)
SODPLOT4:=$(SODSHOCKDIR)$(SODFILE4)

.PHONY: all 

all: sod_plots

sod_plots: $(SODPLOT1) $(SODPLOT2) $(SODPLOT3) $(SODPLOT4)

$(SODPLOT1): $(TESTSDIR)sod_shock_relativistic_problem_1.py
	python $(TESTSDIR)sod_shock_relativistic_problem_1.py

$(SODPLOT2): $(TESTSDIR)sod_shock_relativistic_problem_2.py
	python $(TESTSDIR)sod_shock_relativistic_problem_2.py

$(SODPLOT3): $(TESTSDIR)sod_shock_relativistic_problem_3.py
	python $(TESTSDIR)sod_shock_relativistic_problem_3.py

$(SODPLOT4): $(TESTSDIR)sod_shock_nonrelativistic_problem.py
	python $(TESTSDIR)sod_shock_nonrelativistic_problem.py
