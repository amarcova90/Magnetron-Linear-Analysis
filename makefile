CXX := g++
CXXFLAGS  := -std=c++11 -O3 -W -Wall -Wextra #-Wfatal-errors

TARGET := main
OBJS := main.o nonlinearmodel.o
INCS := util.hpp linearmodel.hpp nonlinearmodel.hpp 
INCLUDES += -I /mnt/c/Users/Andrea/Google\ Drive/Stanford/Research/Magnetron/Magnetron-Linear-Analysis/MTL-4.0.9555-Linux/usr/include

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

%.o: %.cpp $(INCS)
	$(CXX)  $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(OBJS) $(TARGET) *~ *.o solution.txt \
  w_solution.txt solution_figures/*.eps solution_figures/*.png \
  solution_movies/*mp4