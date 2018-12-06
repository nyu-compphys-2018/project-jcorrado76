
PLOTSDIR=../plots/
SODSHOCKDIR=sod_shock_tube_plots/
# the output of savefig puts those single quotes on the filenames 
# need to escape spaces in order to avoid resursive make calling itself 
SODFILE1=\'Low\ Order\ in\ Space\ Low\ order\ in\ Time\ Relativistic\ Sod\ Problem\ 1.jpg\'
SODFILE2=\'Low\ Order\ in\ Space\ Low\ order\ in\ Time\ Relativistic\ Sod\ Problem\ 2.jpg\'
SODFILE3=\'Low\ Order\ in\ Space\ Low\ order\ in\ Time\ Relativistic\ Sod\ Problem\ 3.jpg\'

# need to use := for simply scanned variable expansion equivalent to cocnatenation
SODPLOT1:=$(PLOTSDIR)$(SODSHOCKDIR)$(SODFILE1)
SODPLOT2:=$(PLOTSDIR)$(SODSHOCKDIR)$(SODFILE2)
SODPLOT3:=$(PLOTSDIR)$(SODSHOCKDIR)$(SODFILE3)

# MAKE ALL PLOTS












.PHONY: all sod_plots

all: sod_plots

sod_plots: $(SODPLOT1) $(SODPLOT2) $(SODPLOT3)
	@echo "Regenerating sod plots..."

$(SODPLOT1): sod_shock_relativistic_problem_1.py
	python sod_shock_relativistic_problem_1.py

$(SODPLOT2): sod_shock_relativistic_problem_2.py
	python sod_shock_relativistic_problem_2.py

$(SODPLOT3): sod_shock_relativistic_problem_3.py
	python sod_shock_relativistic_problem_3.py
