CXX := g++
CXXFLAGS := -Wall -Wconversion -Wextra -Wpedantic -fsanitize=address -std=c++11


main: main.cpp
	$(CXX) $(CXXFLAGS) -o main main.cpp

.PHONY: clean
clean:
	$(RM) main
