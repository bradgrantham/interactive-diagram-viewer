CXXFLAGS        +=      -Wall -O

all: lineart graph2lineart

lineart: lineart.cpp
	g++ -L/opt/local/lib  -lglfw -framework OpenGL -framework Cocoa -framework IOkit -lfreeimageplus -lglut lineart.cpp -o lineart -I/opt/local/include
