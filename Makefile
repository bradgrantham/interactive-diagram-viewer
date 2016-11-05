CXXFLAGS        +=      --std=c++11 -Wall -O

all: lineart graph2lineart

lineart: lineart.cpp
	g++ $(CXXFLAGS) -L/opt/local/lib  -lglfw -framework OpenGL -framework Cocoa -framework IOkit -lfreeimageplus lineart.cpp -o lineart -I/opt/local/include
