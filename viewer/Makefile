# Build the MS/MS spectral viewer.
#
#   make            extract spectra, then build the viewer
#   make spectra    extraction only (slow; re-reads every mzML)
#   make viewer     rebuild the HTML from an existing spectra.json
#   make clean      remove generated files
#
# Requires: data/mzml/, data/peak_boundaries.csv, data/Supplemental_Table_2.xlsx

DATA   := ../data
SPECTRA := $(DATA)/spectra.json
OUTPUT  := viewer_MSMS.html

.PHONY: all spectra viewer clean
all: viewer

spectra $(SPECTRA): $(DATA)/peak_boundaries.csv
	python extract_spectra.py \
	    --peaks    $(DATA)/peak_boundaries.csv \
	    --mzml-dir $(DATA)/mzml \
	    --output   $(SPECTRA)

viewer $(OUTPUT): $(SPECTRA) $(DATA)/Supplemental_Table_2.xlsx viewer_template.html
	python build_viewer.py \
	    --spectra  $(SPECTRA) \
	    --table    $(DATA)/Supplemental_Table_2.xlsx \
	    --template viewer_template.html \
	    --output   $(OUTPUT)

clean:
	rm -f $(SPECTRA) $(OUTPUT)
