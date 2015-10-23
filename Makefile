ifndef PVFMM_DIR
$(error Cannot find file: PVFMM_DIR)
endif

ifndef TBSLAS_DIR
$(error Cannot find file: TBSLAS_DIR)
endif

-include $(PVFMM_DIR)/MakeVariables

PSC_INC = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include 
PSC_LIB = -L$(PETSC_DIR)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

TBSLAS_INC = -I$(TBSLAS_DIR)/src -I$(TBSLAS_DIR)/examples/include
CXXFLAGS_TBSLAS = -std=c++0x 

RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include

TARGET_BIN = \
       $(BINDIR)/stokes

all : $(TARGET_BIN)

ifeq ($(INTEL_OFFLOAD_OK),yes) # Have MIC

$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS_PVFMM) $(CXXFLAGS_TBSLAS) -D__MIC_ASYNCH__ $^       $(TBSLAS_INC) $(PSC_LIB) $(LDFLAGS_PVFMM) -o $@
	$(CXX_PVFMM) $(CXXFLAGS_PVFMM) $(CXXFLAGS_TBSLAS) -no-offload      $^_nomic $(TBSLAS_INC) $(PSC_LIB) $(LDFLAGS_PVFMM) -o $@_nomic

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS_PVFMM) $(CXXFLAGS_TBSLAS) -D__MIC_ASYNCH__ $(TBSLAS_INC) $(PSC_INC) -I$(INCDIR) -c $^ -o $@
	$(CXX_PVFMM) $(CXXFLAGS_PVFMM) $(CXXFLAGS_TBSLAS) -no-offload      $(TBSLAS_INC) $(PSC_INC) -I$(INCDIR) -c $^ -o $@_nomic

else
ifeq ($(NVCC_PVFMM),no) # No MIC, No GPU

$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS_PVFMM) $(CXXFLAGS_TBSLAS)                  $^   $(PSC_LIB) $(LDFLAGS_PVFMM) -o $@

else # Have GPU

$(BINDIR)/%: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS_PVFMM) $(CXXFLAGS_TBSLAS)                  $^   $(PSC_LIB) $(LDFLAGS_PVFMM) -o $@
	mv $@ $@_gpu
	touch $@
	chmod +x $@

endif

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX_PVFMM) $(CXXFLAGS_PVFMM) $(CXXFLAGS_TBSLAS)                  $(TBSLAS_INC) $(PSC_INC) -I$(INCDIR) -c $^ -o $@

endif

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) *~ */*~

