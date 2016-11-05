#ifndef __GUI_GLTEXT_H__
#define __GUI_GLTEXT_H__

#include <map>
#include <iostream>
#include <algorithm>

#include <ft2build.h>
#include <freetype/ftbitmap.h>
#include FT_FREETYPE_H

#include <GL/glew.h>

#include "vectormath.h"

namespace GUI {

typedef int fontunits;

inline int topi(int fu) { return fu / 64; }
inline int tofu(int pi) { return pi * 64; }

struct GlyphInfo
{
    int _valid;
    int _index;
    int _texX;
    int _texY;
    int _texW;
    int _texH;
    fontunits _advance;
    int _offsetX;
    int _offsetY;
    GlyphInfo() :
        _valid(false)
    {}
};

struct GLTextureInfo
{
    GLuint id;
    int width;
    int height;
};

struct Face
{
    FT_Face _face;
    int _flags;
    GlyphInfo _charGlyphs[256];
    fontunits _size;
    fontunits _height;
    int _tabWidth;
    GLTextureInfo _textureInfo;
    const GLTextureInfo& GetTextureInfo() { return _textureInfo; }
    void PopulateGlyphs();
    void SizeGlyphMap(int& minw, int& minh);
    void PaintGlyphMap(unsigned char *map, int pitch, int height);
    int DrawCharacter(int x, int y, int ch);
    int CalcCharacterCoordinates(int x, int y, int ch, vec2f *v, vec2f *t);
    void GetStringSize(const char *txt, int *w, int *h);
public:
    void ApplyTexture();
    Face(std::string filename, fontunits units, int flags);
};

// Face manager
struct FaceManager
{
    struct FaceRec
    {
        std::string _filename;
        fontunits _units;
        int _flags;
        FaceRec(std::string filename, fontunits units, int flags) :
            _filename(filename),
            _units(units),
            _flags(flags)
        {
        }
        bool operator<(const FaceRec& f) const
        {
            return (_filename < f._filename) || (_units < f._units) || (_flags < f._flags);
        }
    };
    static std::map<FaceRec, Face*> *_faceMap;
    static Face *GetFace(const std::string& filename, fontunits points, int flags);
};

struct GLText
{
    int _width, _height;
    std::string _str;
    size_t _length;
    Face *_face;
    vec2f *_texcoords;
    vec2f *_vertices;
    GLText(Face *face, const std::string &str);
    void SetText(const std::string& str);
    void Draw(float now); // "float now" so this can subclass Drawable
    int GetWidth() const { return _width; }
    int GetHeight() const { return _height; }
};

} // namespace GUI

using namespace GUI;


#endif /*__GUI_GLTEXT_H__ */
