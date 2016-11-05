#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>

#include <inttypes.h>

#include <map>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>

#include "gltext.h"

namespace GUI {

using namespace std;
using namespace boost;

FT_Library library;

bool gUseKerning = true;

static void InitializeFreeType() __attribute__((constructor));
static void InitializeFreeType()
{
    FT_Init_FreeType(&library);
}

void Face::ApplyTexture()
{
    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();
    glScalef(1.0f / _textureInfo.width, 1.0f / _textureInfo.height, 1);

    glBindTexture(GL_TEXTURE_2D, _textureInfo.id);
}

void Face::PopulateGlyphs()
{
    for(int i = 0; i < 256; i++) {
        int index = FT_Get_Char_Index(_face, i);
        // cout << format("char %d yields glyph %d\n") % i % index;
        if(index != 0) {
            FT_Load_Glyph(_face, index, 0);
            if(_face->glyph->format != FT_GLYPH_FORMAT_BITMAP) {
                FT_Render_Glyph(_face->glyph, FT_RENDER_MODE_NORMAL);
            }

            GlyphInfo &gi = _charGlyphs[i];
            gi._index = index;
            gi._valid = true;
            gi._texW = _face->glyph->bitmap.width;
            gi._texH = _face->glyph->bitmap.rows;
            gi._offsetX = _face->glyph->bitmap_left;
            gi._offsetY = _face->glyph->bitmap_top;
            gi._advance = _face->glyph->advance.x;
        }
    }
}

void Face::SizeGlyphMap(int& minw, int& minh)
{
    int packWidth = 64;
    int packHeight = 1e9;

    do {
        packWidth += 64;
        int curX = 0, curY = 0;
        int rowHeight = 0;

        for(int i = 0; i < 256; i++) {
            GlyphInfo &gi = _charGlyphs[i];
            if(!gi._valid)
                continue;

            // cout << format("%d -> %d\n") % i % gi._texW;

            if((curX + gi._texW) > packWidth) {
                curY += rowHeight + 1;
                curX = 0;
                rowHeight = 0;
            }

            gi._texX = curX;
            gi._texY = curY;

            curX += gi._texW + 1;
            rowHeight = max(rowHeight, gi._texH);
        }

        packHeight = curY + rowHeight;
    } while(packHeight > packWidth);
    // cout << format("for width %d, height is %d\n") % packWidth % packHeight;
    minw = packWidth;
    minh = packHeight;
}

void Face::PaintGlyphMap(unsigned char *map, int pitch, int height)
{
    for(int i = 0; i < 256; i++) {
        GlyphInfo &gi = _charGlyphs[i];
        if(!gi._valid)
            continue;

        int index = FT_Get_Char_Index(_face, i);
        FT_Load_Glyph(_face, index, 0);
        if(_face->glyph->format != FT_GLYPH_FORMAT_BITMAP) {
            FT_Render_Glyph(_face->glyph, FT_RENDER_MODE_NORMAL);
        }

        for(int j = 0; j < gi._texH; j++) {
            for(int i = 0; i < gi._texW; i++){
                int gpitch = _face->glyph->bitmap.pitch;
                int value = _face->glyph->bitmap.buffer[gpitch * j + i];
                int x = gi._texX + i;
                int y = gi._texY + j;
                assert(x < pitch);
                assert(y < height);
                map[pitch * y + x] = value;
            }
        }
    }
}

Face::Face(string filename, fontunits units, int flags) :
    _flags(flags),
    _size(units)
{
    int result = FT_New_Face(library, filename.c_str(), 0, &_face);

    if(result != 0) {
        fprintf(stderr, "freetype error %d creating \"%s\" at %d units\n", result, filename.c_str(), units);
        exit(EXIT_FAILURE);
    }

    FT_Set_Char_Size(_face, 0, _size, 72, 72);

    _tabWidth = 8 * topi(units) / 2;
    _height = _face->size->metrics.height;

    GLTextureInfo &tex = _textureInfo;

    PopulateGlyphs();
    SizeGlyphMap(tex.width, tex.height);

    unsigned char *map = new unsigned char[tex.width * tex.height];
    for(int i = 0; i < tex.width * tex.height; map[i] = 0, i++);

    PaintGlyphMap(map, tex.width, tex.height);

    glGenTextures(1, &tex.id);
    glBindTexture(GL_TEXTURE_2D, tex.id);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA8, tex.width, tex.height, 0, GL_ALPHA, GL_UNSIGNED_BYTE, map);

    glGenerateMipmapEXT(GL_TEXTURE_2D);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glBindTexture(GL_TEXTURE_2D, 0);

    delete[] map;
}

int Face::DrawCharacter(int x, int y, int ch)
{
    GlyphInfo &gi = _charGlyphs[ch];

    if(!gi._valid)
        return 0;

    int w = gi._texW;
    int h = gi._texH;
    int tx = gi._texX;
    int ty = gi._texY;
    int dx = gi._offsetX;
    int dy = gi._offsetY - h;

    glBegin(GL_QUADS);
        glTexCoord2f(tx, ty + h);
        glVertex2f(x + dx, y + dy);

        glTexCoord2f(tx + w, ty + h);
        glVertex2f(x + dx + w, y + dy);

        glTexCoord2f(tx + w, ty);
        glVertex2f(x + dx + w, y + dy + h);

        glTexCoord2f(tx, ty);
        glVertex2f(x + dx, y + dy + h);
    glEnd();

    return gi._advance / 64;
}

int Face::CalcCharacterCoordinates(int x, int y, int ch, vec2f *v, vec2f *t)
{
    GlyphInfo &gi = _charGlyphs[ch];

    if(!gi._valid)
        return 0;

    int w = gi._texW;
    int h = gi._texH;
    int tx = gi._texX;
    int ty = gi._texY;
    int dx = gi._offsetX;
    int dy = gi._offsetY - h;

    t[0][0] = tx;
    t[0][1] = ty + h;
    v[0][0] = x + dx;
    v[0][1] = y + dy;

    t[1][0] = tx + w;
    t[1][1] = ty + h;
    v[1][0] = x + dx + w;
    v[1][1] = y + dy;

    t[2][0] = tx + w;
    t[2][1] = ty;
    v[2][0] = x + dx + w;
    v[2][1] = y + dy + h;

    t[3][0] = tx;
    t[3][1] = ty;
    v[3][0] = x + dx;
    v[3][1] = y + dy + h;

    return gi._advance / 64;
}

void Face::GetStringSize(const char *txt, int *w, int *h)
{
    int posx = 0;
    const char *p;
	int localx = posx;
	bool kerning_valid = false;
	int prev_char;
    FT_Face face = _face;
    int miny = 100000, maxy = -100000;

	for(p = txt; *p; p++) {

		if(*p == '\t') {
			int tab_pos = localx + tofu(_tabWidth) - posx;
			tab_pos = tab_pos / tofu(_tabWidth) * tofu(_tabWidth);
			localx = posx + tab_pos;
			kerning_valid = false;
			continue;
		}

        GlyphInfo &gi = _charGlyphs[(int)*p];

		if(!gi._valid) {
			kerning_valid = false;
			continue;
		}

		if(gUseKerning && kerning_valid) {
			FT_Vector delta;
			FT_Get_Kerning(face, prev_char, gi._index, 0, &delta);
			if(delta.x != 0)
                localx += delta.x;
		}

        localx += gi._advance;
        miny = min(miny, gi._offsetY);
        maxy = max(maxy, gi._offsetY + gi._texH);

		kerning_valid = true;
		prev_char = gi._index;
	}
    *w = topi(localx);
    *h = maxy - miny; // topi(face->size->metrics.height);
}

map<FaceManager::FaceRec, Face*> *FaceManager::_faceMap;

Face *FaceManager::GetFace(const string& filename, fontunits points, int flags)
{
    map<FaceRec, Face*>::iterator i;
    Face *face;
    if(_faceMap == NULL)
        _faceMap = new map<FaceRec, Face*>();

    i = _faceMap->find(FaceRec(filename, points, flags));
    if(i == _faceMap->end()) {
        face = new Face(filename, points, flags);
        (*_faceMap)[FaceRec(filename, points, flags)] = face;
    } else {
        face = (*_faceMap)[FaceRec(filename, points, flags)];
    }
    return face;
}

GLText::GLText(Face *face, const string &str) :
    _face(face),
    _texcoords(NULL),
    _vertices(NULL)
{
    SetText(str);
}

void GLText::SetText(const string &str)
{
    delete _texcoords;
    delete _vertices;
    _str = str;
    _length = _str.length();
    _texcoords = new vec2f[4 * _length];
    _vertices = new vec2f[4 * _length];

    _face->GetStringSize(str.c_str(), &_width, &_height);

	int localx = 0;
	bool kerning_valid = false;
	int prev_char;
    FT_Face face = _face->_face;

    int i = 0;
    _str = str;
	for(const char *p = str.c_str(); *p; i++, p++) {

		if(*p == '\t') {
			int tab_pos = localx + tofu(_face->_tabWidth);
			tab_pos = tab_pos / tofu(_face->_tabWidth) * tofu(_face->_tabWidth);
			localx = tab_pos;
			kerning_valid = false;
			continue;
		}

        GlyphInfo &gi = _face->_charGlyphs[(int)*p];

		if(!gi._valid) {
			kerning_valid = false;
			continue;
		}

		if(gUseKerning && kerning_valid) {
			FT_Vector delta;
			FT_Get_Kerning(face, prev_char, gi._index, 0, &delta);
			if(delta.x != 0)
                localx += delta.x;
		}

        _face->CalcCharacterCoordinates(topi(localx), 0, *p, &_vertices[4 * i], &_texcoords[4 * i]);

        localx += gi._advance;

		kerning_valid = true;
		prev_char = gi._index;
	}
}

void GLText::Draw(float now)
{
    glEnable(GL_TEXTURE_2D);
    _face->ApplyTexture();
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, _vertices);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glTexCoordPointer(2, GL_FLOAT, 0, _texcoords);
    glDrawArrays(GL_QUADS, 0, _length * 4);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_2D);
}

} // namespace GUI

// vi: ts=4
