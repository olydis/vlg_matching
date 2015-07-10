include ../Make.helper
CXX_FLAGS = $(MY_CXX_FLAGS) # in compile_options.config
LIBS = -lsdsl 
SRC_DIR = src
TMP_DIR = ../tmp

# TODO: configure those:
SAMPLES = $(shell seq 1000 1000 20000)
COLL_ID = sources

TC_IDS:=$(call config_ids,test_case.config)

all: execs 

GM_WORDS = collections/$(COLL_ID).words.txt
GM_PATTERNS = $(foreach SAMPLE,$(SAMPLES),collections/$(COLL_ID).$(SAMPLE).txt)
GM_EXECS = $(foreach TC_ID,$(TC_IDS),build/gm_index-$(TC_ID).x) \
           $(foreach TC_ID,$(TC_IDS),build/gm_search-$(TC_ID).x)
INDEX_SENTINELS = $(foreach TC_ID,$(TC_IDS),collections/$(COLL_ID).$(TC_ID).sentinel)

RES_FILES = $(foreach SAMPLE,$(SAMPLES),$(foreach TC_ID,$(TC_IDS),results/$(TC_ID).$(SAMPLE)))
						
RES_FILE=results/all.txt

# Targets for index and search
# TODO

execs: $(GM_EXECS)

timing: execs $(RES_FILES) 
	cat $(RES_FILES) > $(RES_FILE)
	@cd visualize; make	

# Format: results/[TC_ID].[SAMPLE]
results/%: $(GM_PATTERNS)
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval SAMPLE:=$(call dim,2,$*))
	@echo "# COLL_ID = $(COLL_ID)" > $@
	@echo "# PATT_SAMPLE = $(SAMPLE)" >> $@
	@echo "# TC_ID = $(TC_ID)" >> $@
	$(eval SENTINEL:="collections/$(COLL_ID).$(call dim,1,$*).sentinel")
	@test -e $(SENTINEL) || { echo "Running build/gm_index-$(TC_ID).x -c collections/$(COLL_ID)"; ./build/gm_index-$(TC_ID).x -c collections/$(COLL_ID); touch $(SENTINEL); }
	@echo "Running build/gm_search-$(TC_ID).x -c collections/$(COLL_ID) -p collections/$(COLL_ID).$(SAMPLE).txt"
	@build/gm_search-$(TC_ID).x -c collections/$(COLL_ID) -p collections/$(COLL_ID).$(SAMPLE).txt >> $@

$(GM_WORDS):
	@echo "Generating words"
	@../../examples/generate-pattern.x collections/$(COLL_ID)/$(COLL_ID).raw 20 4 > $(GM_WORDS)

collections/$(COLL_ID).%.txt: $(GM_WORDS)
	$(eval SAMPLE:=$(call dim,1,$*))
#	$(eval MINGAP:=$(SAMPLE))
#	$(eval MAXGAP:=$(shell echo $(MINGAP)+10 | bc))
	$(eval MINGAP:=0)
	$(eval MAXGAP:=$(SAMPLE))
	@sed 'N;s/\n/.{$(MINGAP),$(MAXGAP)}/' $(GM_WORDS) > $@

# TODO: less hacky
build/%:
	cd build && cmake .. && make -j

clean-build:
	@echo "Remove executables"
	rm -f $(GM_EXECS)

clean:
	rm -f $(GM_EXECS)
	rm -f $(GM_PATTERNS)
	rm -f $(GM_WORDS)
	rm -f $(INDEX_SENTINELS)

clean_results: 
	rm -f $(RES_FILES) 
		

cleanall: clean clean_results
