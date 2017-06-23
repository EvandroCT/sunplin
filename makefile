#in order to build device code, nvidia compiler must be invoked
CC=nvcc

#headers location
INCDIR=include/

#source location
SRCDIR=src/

#R package source location
PKGDIR=$(SRCDIR)rpkg/

#builded libraries location
LIBDIR=lib/

#required R includes location
RINC=$(shell R CMD config --cppflags)

#source extension
SRCEXT=cu

#standalone app source files
STDSOURCES=$(shell ls $(SRCDIR)*.$(SRCEXT))
#$(info $$var is [${STDSOURCES}])
#R package source files
PKGSOURCES=$(shell ls $(PKGDIR)*.$(SRCEXT))

#builded objects location
BUILDIR=build/

#list of builded objects
OBJS=$(patsubst $(SRCDIR)%,$(BUILDIR)%,$(STDSOURCES:.$(SRCEXT)=.o))

#list of headers
INCLUDES=$(shell ls $(INCDIR)*)	

#standalone app target
STDTARGET=sunplin

#R package target
PKGTARGET=$(LIBDIR)librsunplin.so

#compiler flags (dc flag allows device code to be linked)
CFLAGS=-arch=sm_35 -std=c++11 -Xcompiler -fPIC -w -I $(INCDIR) -dc

#flags to generate R's loadable shared object (rdc flag allows device code to be linked)
SOFLAGS=-arch=sm_35 -std=c++11 -Xcompiler -fPIC -w -I $(INCDIR) -shared $(RINC) -rdc=true

#linker flags
LDFLAGS=-arch=sm_35

all: $(STDTARGET) $(PKGTARGET)

rpackage: rpkg

rpkg: $(PKGTARGET)
	
$(STDTARGET): $(OBJS)
	@echo "Linking..."
	$(CC) $(LDFLAGS) $(OBJS) -o $(STDTARGET)

$(PKGTARGET): $(OBJS) $(PKGSOURCES) $(INCLUDES)
	@echo "Creating shared object..."
	@mkdir -p $(LIBDIR);
	$(CC) $(SOFLAGS) $(PKGSOURCES) $(OBJS) -o $(PKGTARGET)

$(BUILDIR)%.o: $(SRCDIR)%.$(SRCEXT) $(INCLUDES)
	@echo "Compiling..."
	@mkdir -p $(BUILDIR);
	$(CC) $(CFLAGS) $< -o $@

clean:
	@echo "Cleaning... "
	@rm -Rf $(BUILDIR) $(LIBDIR) sunplin
