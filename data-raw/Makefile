UNAME ?= $(shell uname)
ifeq ($(UNAME),Darwin)
	IJ ?= /Applications/Fiji.app/Contents/MacOS/ImageJ-macosx
	JAVA_HOME ?= /Library/Java/JavaVirtualMachines/openjdk-11.jdk/Contents/Home
else ifeq ($(UNAME),Linux)
	IJ ?= ImageJ-linux64
	JAVA_HOME ?= $(shell $(IJ) --print-java-home)
endif
IJ_ARGS ?= --java-home $(JAVA_HOME) --ij2 --headless
IJM ?= extract_data_points_from_image_plot-single.ijm

INPUT_ROISET ?= \
nejm199501263320402_f1-subject1-RoiSet.zip \
nejm199501263320402_f1-subject2-RoiSet.zip \
nejm199501263320402_f1-subject3-RoiSet.zip \
nejm199501263320402_f1-subject4-RoiSet.zip \
nejm199501263320402_f1-subject5-RoiSet.zip \
nejm199501263320402_f1-subject6-RoiSet.zip
OUTPUT_CSV := $(INPUT_ROISET:zip=csv)

ROOT_DIR := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))


../data/pantaleo1995-figure1.csv : $(OUTPUT_CSV)
	awk -F, 'BEGIN {print "group,year,cd4_cells_per_mm3"} FNR>1 {print substr(FILENAME,24,8) "," sprintf("%2.1f",$$1) "," sprintf("%.0f",$$2)}' $^ > $@

$(OUTPUT_CSV) : %.csv : %.zip $(IJM)
	$(IJ) $(IJ_ARGS) --run $(IJM) 'file_input_roiset="$(ROOT_DIR)/$<"'
