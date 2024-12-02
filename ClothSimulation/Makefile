#Makefile for cloth simulation
ASSIGN    = cloth_simulation
BREWPATH  = $(shell brew --prefix)
CXX       = /opt/homebrew/opt/llvm/bin/clang++
CXXFLAGS  = -std=c++11 -D_CRT_SECURE_NO_WARNINGS -DGL_SILENCE_DEPRECATION -Wno-macro-redefined -O2 -I./Headers -I$(BREWPATH)/include -I./imgui -I./imgui/backends -Xpreprocessor -fopenmp -arch arm64
LDFLAGS   = -L$(BREWPATH)/lib -lglfw -framework OpenGL -framework Cocoa -framework IOKit -L/opt/homebrew/opt/libomp/lib -lomp -arch arm64

# ImGui 源文件
IMGUI_SRCS = imgui/imgui.cpp imgui/imgui_draw.cpp imgui/imgui_tables.cpp imgui/imgui_widgets.cpp \
             imgui/backends/imgui_impl_glfw.cpp imgui/backends/imgui_impl_opengl3.cpp
IMGUI_OBJS = $(IMGUI_SRCS:.cpp=.o)

all: $(ASSIGN)

$(ASSIGN): main.o glad.o $(IMGUI_OBJS)
	$(CXX) $^ $(LDFLAGS) -o $@

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

glad.o: glad.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(ASSIGN) *.o imgui/*.o imgui/backends/*.o *~ *.dSYM

.PHONY: clean all