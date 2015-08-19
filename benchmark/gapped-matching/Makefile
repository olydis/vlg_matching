include ../Make.helper
CXX_FLAGS = $(MY_CXX_FLAGS) # in compile_options.config
LIBS = -lsdsl 
SRC_DIR = src
TMP_DIR = ../tmp

.SECONDARY:

TC_IDS:=$(call config_ids,test_case.config)
UNIQUE_SP_LENGTHS:=$(shell echo $(call config_column,test_case.config,2) | tr " " "\n" | uniq | tr "\n" " ")

COLL_ID = sources-1048576

ALGOS:=$(call config_ids,algorithms.config)

all: execs 

GM_WORDS = $(foreach SP_LENGTH,$(UNIQUE_SP_LENGTHS),collections/$(COLL_ID).$(SP_LENGTH).words.txt)
GM_PATTERNS = $(foreach TC_ID,$(TC_IDS),collections/$(COLL_ID).$(TC_ID).regex.txt)
GM_EXECS = $(foreach ALGO,$(ALGOS),build/gm_index-$(ALGO).x) \
           $(foreach ALGO,$(ALGOS),build/gm_search-$(ALGO).x)
INDEX_SENTINELS = $(foreach ALGO,$(ALGOS),collections/$(COLL_ID).$(ALGO).sentinel)

RES_FILES = $(foreach TC_ID,$(TC_IDS),$(foreach ALGO,$(ALGOS),results/$(ALGO).$(TC_ID)))
						
RES_FILE=results/all.txt

execs: $(GM_EXECS)

timing: execs $(RES_FILES) 
	cat $(RES_FILES) > $(RES_FILE)
	@cd visualize; make	

# Format: results/[ALGO].[TC_ID]
results/%: $(GM_PATTERNS)
	@mkdir -p results
	$(eval ALGO:=$(call dim,1,$*))
	$(eval TC_ID:=$(call dim,2,$*))
	@echo "# COLL_ID = $(COLL_ID)" > $@
	@echo "# PATT_SAMPLE = $(TC_ID)" >> $@
	@echo "# ALGO = $(ALGO)" >> $@
	$(eval SENTINEL:="collections/$(COLL_ID).$(call dim,1,$*).sentinel")
	@test -e $(SENTINEL) || { echo "Running build/gm_index-$(ALGO).x -c collections/$(COLL_ID)"; ./build/gm_index-$(ALGO).x -c collections/$(COLL_ID); touch $(SENTINEL); }
	@echo "Running build/gm_search-$(ALGO).x -c collections/$(COLL_ID) -p collections/$(COLL_ID).$(TC_ID).txt"
	@build/gm_search-$(ALGO).x -c collections/$(COLL_ID) -p collections/$(COLL_ID).$(TC_ID).regex.txt >> $@

collections/%.words.txt:
	$(eval COLL_ID:=$(call dim,1,$*))
	$(eval SP_LENGTH:=$(call dim,2,$*))
	@echo "Extracting words of length $(SP_LENGTH) from $(COLL_ID)..."
	@../../examples/generate-pattern.x collections/$(COLL_ID)/$(COLL_ID).raw 20 $(SP_LENGTH) > $@

collections/%.regex.txt: $(GM_WORDS)
	$(eval COLL_ID:=$(call dim,1,$*))
	$(eval TC_ID:=$(call dim,2,$*))
	$(eval SP_LENGTH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval GAP:=$(call config_select,test_case.config,$(TC_ID),3))
	@echo "$(TC_ID)"
	@echo "Creating regex patterns with words of length $(SP_LENGTH) from $(COLL_ID) and a gap of $(GAP)..."
	@sed 'N;s/\n/.{$(GAP)}/' collections/$(COLL_ID).$(SP_LENGTH).words.txt > $@

# TODO: less hacky
build/%:
	cd build && cmake .. && make -j

clean-build:
	@echo "Remove executables"
	rm -rf build/*

clean:
	rm -f $(GM_EXECS)
	rm -f $(GM_PATTERNS)
	rm -f $(GM_WORDS)
	rm -f $(INDEX_SENTINELS)

clean-results: 
	rm -f $(RES_FILES) 
		

cleanall: clean clean-results clean-build
