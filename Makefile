CXXFLAGS        +=      --std=c++11 -Wall -g -O

all: lineart graph2lineart 

lineart: lineart.cpp gltext.cpp vectormath.cpp
	g++ $(CXXFLAGS) -L/opt/local/lib  -lglfw -lglew -framework OpenGL -framework Cocoa -framework IOkit -lfreeimageplus -lfreetype lineart.cpp gltext.cpp vectormath.cpp -o lineart -I/opt/local/include -I/opt/local/include/freetype2
