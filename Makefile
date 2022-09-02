CC = g++

CFLAGS = -Wall -O3

# define any directories containing header files other than /usr/include
INCLUDES = -I./

# define library paths in addition to /usr/lib
LFLAGS = -L./ -L/usr/local/lib

# define any libraries to link into executable:
LIBS = -lIL -lboost_program_options
# DevIL = IL, CDT does not need to be linked, only source files

# define the C source files
SRCS = main.cpp Document.cpp

OBJS = $(SRCS:.c=.o)

# define the executable file
MAIN = dbrt


.PHONY: depend clean

all:    $(MAIN)
	@echo  Document Background Removing Tool has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)


.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
