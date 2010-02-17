export MAKEFILES_DIR=makefiles

include $(MAKEFILES_DIR)/Makefile.dirs.vars

.PHONY: init depend clean clean_obj libs apps all compile

all: libs apps test


depend: 
	$(MAKE) -I $(MAKEFILES_DIR) -f $(MAKEFILES_DIR)/Makefile.depend depend

libs: depend
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.build libs

apps: depend
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.build apps

compile: depend
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.build compile

test: depend
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.tests test

performance_test: depend
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.tests performance_test

coverage: depend
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.lcov coverage

cppcheck:
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.cppcheck cppcheck

cccc:
	$(MAKE) -I $(MAKEFILES_DIR)  -f $(MAKEFILES_DIR)/Makefile.cccc cccc

clean:
	$(RM) -r $(BUILD_DIR)/*
	$(RM) $(BUILD_DIR)/.init*



# DO NOT DELETE
