#SRC_DIR := .../src
#OBJ_DIR := .../obj
SRC_DIR := src
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
#LDFLAGS := ...
#CPPFLAGS := ...
#CXXFLAGS := ...

#main.exe: $(OBJ_FILES)
	#g++ $(LDFLAGS) -o $@ $^

#$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	#g++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

main.exe: $(OBJ_FILES)
	g++ -g -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ -g -c -o $@ $<

plots: plots/sod_shock_tube.jpg
	python plot_final_state.py

clean: 
	rm $(OBJ_DIR)/*.o 
	rm main.exe
