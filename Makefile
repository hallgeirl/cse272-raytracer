NAME = miro


all: $(NAME)

include Makedefs
SOURCES -=  parse.cpp lexer.cpp
OBJS -=  parse.o lexer.o
PDF=$(shell ls *.pdf 2>/dev/null)
PNG=$(shell ls *.png 2>/dev/null)
PLOTS=$(patsubst %.p,%.pdf,$(shell ls *.p))
#PLOTS=$(patsubst %.p,%.png,$(shell ls *.p))

CXXFLAGS += -DLINUX -DALTERNATIVE -g

assignment2.cpp: assignment2.h

# 
# lexer.cpp: lexer.lex
# 	$(ECHO) "Flex-ing lexer.lex"
# 	$(FLEX) -o$@ lexer.lex
# 
# parse.cpp: parse.y
# 	$(ECHO) "Bison-ing parse.y"
# 	$(BISON) -d -o $@ parse.y
# 	@if [ -f parse.hpp ]; then \
# 		mv parse.hpp parse.cpp.h; \
# 	fi
# 	@if [ -f parse.tab.hpp ]; then \
# 		mv parse.tab.hpp parse.cpp.h; \
# 	fi
# 	@if [ -f parse.tab.h ]; then \
# 		mv parse.tab.h parse.cpp.h; \
# 	fi

-include .deps/*.d

$(NAME): $(OBJS)
	$(ECHO) "Linking $@..."
	$(CXX) -o $@ $(OBJS) $(LDFLAGS)
	$(ECHO) "Built $@!"

freeimage:
	$(MAKE) -C lib/src/FreeImage
	mkdir -p lib/lib
	mv lib/src/FreeImage/Dist/libfreeimage.a lib/lib/libfreeimage.a

plots: $(PLOTS)
	cp *.pdf report/a2/plots
	
test: $(NAME)
	./miro && eog *.ppm && rm *.ppm

