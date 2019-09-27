CXX      := -c++
#CXXFLAGS := -pedantic-errors -std=c++11 -Wall -Wextra -Werror
CXXFLAGS := -pedantic-errors -std=c++11 -w
LDFLAGS  := -L/usr/lib -lstdc++ -lm
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := ./
TARGET   := prog-opt
INCLUDE  := -Iinclude/
SRC      :=                      \
   $(wildcard src/*.C)           \

# $(wildcard src/module1/*.C)
# $(wildcard src/module1/*.C)
# $(wildcard src/module1/*.C)

OBJECTS := $(SRC:%.C=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.C
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(APP_DIR)/$(TARGET) $(OBJECTS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
