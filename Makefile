PREFIX=$(HOME)

LIBNAME = PLASMA
LIBRARY = lib$(LIBNAME).so

CXX      = g++
OPT2     = -O3
CXXFLAGS = $(OPT2) -Wall -fPIC -I$(PREFIX)/include
SOFLAGS  = -shared

SRCS = $(wildcard *.C)
EXES = $(SRCS:.C=.exe)

SOURCES = $(wildcard *.cc)
HEADERS = $(SOURCES:.cc=.h)
OBJECTS = $(SOURCES:.cc=.o)
DEPFILE = $(SOURCES:.cc=.d)


all: $(EXES)

# include *.d files, which are makefiles defining dependencies between files
ifeq ($(filter uninstall clean tags, $(MAKECMDGOALS)),)
  -include $(DEPFILE)
endif

# rules to create *.d files
%.d:%.cc
	@echo creating $@
	@set -e; rm -f $@; \
	  $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$ 

# lib$(LIBNAME).so depends on all *.o files.
#  The flag "-shared" is used to create shared libs
#  $@ represents the target, that is, lib$(LIBNAME).so
#  $^ represents all the prerequisites, i.e., all *.o files
$(LIBRARY): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LIBS) $(SOFLAGS) -o $@ $^

# An xxx.o file depends on xxx.cc. It is created with the command:
# 	g++ -c xxx.cc -o xxx.o
# Since this is obvious, "make" automatically does it. 
# There is no need to explicitly write down the rules to do it.

clean:
	$(RM) *.o *.d *.d.* *~ $(LIBRARY)

tags:
	ctags --c-kinds=+p $(HEADERS) $(SOURCES)

install: $(LIBRARY)
	@echo "PREFIX=$(PREFIX)"
	@echo -n "checking if $(PREFIX) exists..."
	@if [ -d $(PREFIX) ]; then \
	  echo "yes."; \
	else \
	  echo "no."; \
	  echo "mkdir $(PREFIX)..."; \
	  mkdir $(PREFIX); \
	fi
	@echo -n "copying $(LIBRARY) to $(PREFIX)/lib..."
	@if [ -d $(PREFIX)/lib ]; then \
	  cp $(LIBRARY) $(PREFIX)/lib; \
	else \
	  mkdir $(PREFIX)/lib; \
	  cp $(LIBRARY) $(PREFIX)/lib; \
	fi; 
	@echo "done."; 
	@echo -n "copying *.h to $(PREFIX)/include/$(LIBNAME)..."
	@if [ -d $(PREFIX)/include ]; then \
	  if [ -d $(PREFIX)/include/$(LIBNAME) ]; then \
	    cp *.h $(PREFIX)/include/$(LIBNAME); \
	  else \
	    mkdir $(PREFIX)/include/$(LIBNAME); \
	    cp *.h $(PREFIX)/include/$(LIBNAME); \
	  fi; \
	else \
	  mkdir $(PREFIX)/include; \
	  mkdir -p $(PREFIX)/include/$(LIBNAME); \
	  cp *.h $(PREFIX)/include/$(LIBNAME); \
	fi
	@echo "done."

uninstall:
	$(RM) $(PREFIX)/lib/$(LIBRARY)
	$(RM) -r $(PREFIX)/include/$(LIBNAME)

$(EXES):%.exe:%.C install
	$(CXX) $< $(CXXFLAGS) -L. -l$(LIBNAME) $(LIBS) -o $@

.PHONY: all tags clean install uninstall
