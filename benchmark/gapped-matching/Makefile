include ../Make.helper
CXX_FLAGS = $(MY_CXX_FLAGS) # in compile_options.config
LIBS = -lsdsl 
SRC_DIR = src
TMP_DIR = ../tmp

# TODO: configure those:
COLL_ID = ulysses
COLL_PATT_ID = ulysses/patterns/test_patterns.txt

TC_IDS:=$(call config_ids,test_case.config)

all: execs 

GM_EXECS = $(foreach TC_ID,$(TC_IDS),build/gm_index-$(TC_ID).x) \
           $(foreach TC_ID,$(TC_IDS),build/gm_search-$(TC_ID).x)

RES_FILES = $(foreach TC_ID,$(TC_IDS),results/$(TC_ID))
						
RES_FILE=results/all.txt

# Targets for index and search
# TODO

execs: $(GM_EXECS) 


timing: execs $(RES_FILES) 
	cat $(RES_FILES) > $(RES_FILE)
	@cd visualize; make	

# Format: results/[TC_ID]
results/%: 
	$(eval TC_ID:=$(call dim,1,$*)) 
	@echo "Running build/gm_index-$(TC_ID).x -c collections/$(COLL_ID)"
	@echo "Running build/gm_search-$(TC_ID).x -c collections/$(COLL_ID) -p collections/$(COLL_PATT_ID)"
	@echo "# COLL_ID = $(COLL_ID)" >> $@
	@echo "# COLL_PATT_ID = $(COLL_PATT_ID)" >> $@
	@echo "# TC_ID = $(TC_ID)" >> $@
	@build/gm_index-$(TC_ID).x -c collections/$(COLL_ID) >> $@
	@build/gm_search-$(TC_ID).x -c collections/$(COLL_ID) -p collections/$(COLL_PATT_ID) >> $@

clean-build:
	@echo "Remove executables"
	rm -f $(GM_EXECS)

clean:
	rm -f $(GM_EXECS)

clean_results: 
	rm -f $(RES_FILES) 
		

cleanall: clean clean_results
