CXX := g++
CXXFLAGS  := -std=c++11 -O3 -W -Wall -Wextra #-Wfatal-errors 

current_dir := $(shell pwd)

OBJS_LIN := main_linear.o
OBJS_VAL := main_NLvalidation.o nonlinearmodel.o
OBJS_NL  := main_nonlinear.o nonlinearmodel.o
INCS_LIN := util.hpp linearmodel.hpp
INCS_VAL := util.hpp linearmodel.hpp nonlinearmodel.hpp 
INCS_NL := util.hpp linearmodel.hpp nonlinearmodel.hpp 
INCLUDES = -I $(current_dir)/MTL-4.0.9555-Linux/usr/include -lboost_system -fopenmp
TARGET_LIN := run_linear
TARGET_VAL := run_NLvalidation
TARGET_NL := run_NL

.PHONY: all linear validation nonlinear

all : linear validation nonlinear
linear : $(TARGET_LIN)
validation : $(TARGET_VAL)
nonlinear : $(TARGET_NL)

$(TARGET_LIN): $(OBJS_LIN)
	$(CXX) $(OBJS_LIN) -o $(TARGET_LIN)

$(TARGET_VAL): $(OBJS_VAL)
	$(CXX) $(OBJS_VAL) $(INCLUDES) -o $(TARGET_VAL)

$(TARGET_NL): $(OBJS_NL)
	$(CXX) $(OBJS_NL) $(INCLUDES) -o $(TARGET_NL)

%.o: %.cpp $(INCS_LIN) $(INCS_VAL) $(INCS_NL)
	$(CXX)  $(CXXFLAGS) $(INCLUDES) -o $@ -c $< 

.PHONY: clean
clean:
	$(RM) $(OBJS_LIN) $(TARGET_NL) $(TARGET_LIN) $(TARGET_VAL) *~ *.o solution.txt \
  w_solution.txt solution_figures/*.eps solution_figures/*.png \
  solution_movies/*mp4 validation_files/*.txt validation_files/*.png \
  linear_model_solutions/*.txt