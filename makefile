CXX := g++
CXXFLAGS := -O3 -Wall -Wconversion -Wextra -Wpedantic -std=c++11 

TARGET := main
OBJS := main.o
INCS := util.hpp linearmodel.hpp

$(TARGET): $(OBJS)
	$(CXX)  -o $(TARGET) $(OBJS)

%.o: %.cpp $(INCS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

.PHONY: clean
clean:
	$(RM) $(OBJS) $(TARGET) *~ *.o solution.txt \
  w_solution.txt solution_figures/*.eps solution_figures/*.png \
  solution_movies/*mp4