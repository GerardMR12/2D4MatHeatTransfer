# Compiler settings
CXX = g++
CXXFLAGS = -Wall -std=c++11

# Directories
SCRIPTS_DIR = scripts
SRC_DIR = src
BUILD_DIR = build

# Target executable name
TARGET = $(BUILD_DIR)/MyApp
PYTHON = python

# Create directory file
$(shell mkdir $(BUILD_DIR))

# Default target
all: $(TARGET) post-process clean
	@echo ----- Program finished!

# Link the target with objects
$(TARGET): $(BUILD_DIR)/main.o
	@echo ----- Linking the target...
	@$(CXX) $(CXXFLAGS) -o $(TARGET) $(BUILD_DIR)/main.o
	@echo ----- Running the target...
	@./$(TARGET)

# Compile the source files into object files
$(BUILD_DIR)/main.o: src/main.cpp
	@echo ----- Compiling the source files...
	@$(CXX) $(CXXFLAGS) -c src/main.cpp -o $(BUILD_DIR)/main.o

# Run Python script for post-processing
post-process:
	@echo ----- Running python plotting script...
	@$(PYTHON) $(SCRIPTS_DIR)/plot.py

# Clean target
clean:
	@echo ----- Cleaning environment...
	@powershell -Command "Remove-Item -Path '$(BUILD_DIR)' -Force -Recurse"
	@if exist $(BUILD_DIR) rmdir /S /Q $(BUILD_DIR)

.PHONY: all post-process clean