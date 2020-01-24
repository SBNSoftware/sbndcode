DEPS        = BranchType.hh MetricHolder.cc NueRecoOptimiser.hh
OBJFILES    = MetricHolder.o NueRecoOptimiser.o MakeSelectionEfficiencyPlots.o
CFILES      = MetricHolder.cc NueRecoOptimiser.cc MakeSelectionEfficiencyPlots.cc
EXECUTABLE  = MakeSelectionEfficiencyPlots

LIBS        = `root-config --glibs`
SOFLAGS     = -fPIC -shared
ROOTCFLAGS:= `root-config --cflags`
ROOTLFLAGS:= `root-config --glibs`
CXX         = `root-config --cxx`
SOFLAGS     = -fPIC -shared
LIBS        = `root-config --glibs`
CXXFLAGS    = `root-config --cflags`

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJFILES)
	@echo Building Exectuable $@
	@$(CXX) $(CXXFLAGS) -Wall -lTMVA -lTMVAGui -o $@ $^ $(LIBS) 

%.o : %.cc $(DEPS)
	@echo Building Object File $@
	@$(CXX) $(CXXFLAGS) -c $^ $(LIBS)

.PRECIOUS: $(OBJFILES)

.PHONY: clean

clean:
	@rm -f ./*~ $(OBJFILES) $(EXECUTABLES)

