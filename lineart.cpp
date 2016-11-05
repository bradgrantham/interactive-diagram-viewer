#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <cstring>
#include <cfloat>
#include <cstdarg>
#include <GL/glut.h>

#include <vector>
#include <string>

using namespace std;


float vec3fOrigin[3] = {0.0f, 0.0f, 0.0f}; 
float vec3fXAxis[3] = {1.0f, 0.0f, 0.0f}; 
float vec3fYAxis[3] = {0.0f, 1.0f, 0.0f}; 
float vec3fZAxis[3] = {0.0f, 0.0f, 1.0f}; 

float vec4fOrigin[4] = {0.0f, 0.0f, 0.0f, 1.0f}; 
float vec4fXAxis[4] = {1.0f, 0.0f, 0.0f, 0.0f}; 
float vec4fYAxis[4] = {0.0f, 1.0f, 0.0f, 0.0f}; 
float vec4fZAxis[4] = {0.0f, 0.0f, 1.0f, 0.0f}; 


void vec3fCopy(float src[3], float dest[3])
{
    dest[0] = src[0];
    dest[1] = src[1];
    dest[2] = src[2];
}

float vec3fDot(float v0[3], float v1[3])
{
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}

void vec3fCross(float v0[3], float v1[3], float result[3])
{
    result[0] = v0[1] * v1[2] - v0[2] * v1[1];
    result[1] = v0[2] * v1[0] - v0[0] * v1[2];
    result[2] = v0[0] * v1[1] - v0[1] * v1[0];
}

void vec3fBlend(float v0[3], float w0, float v1[3], float w1, float result[3])
{
    result[0] = v0[0] * w0 + v1[0] * w1;
    result[1] = v0[1] * w0 + v1[1] * w1;
    result[2] = v0[2] * w0 + v1[2] * w1;
}

void vec3fScale(float v[3], float w, float r[3])
{
    r[0] = v[0] * w;
    r[1] = v[1] * w;
    r[2] = v[2] * w;
}

void vec3fNormalize(float v[3], float r[3])
{
    vec3fScale(v, 1.0 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]), r);
}

void vec3fReflect(float i[3], float n[3], float r[3])
{
    float dot;

    dot = vec3fDot(i, n);

    r[0] = i[0] - 2.0f * n[0] * dot;
    r[1] = i[1] - 2.0f * n[1] * dot;
    r[2] = i[2] - 2.0f * n[2] * dot;
}

void vec4fCopy(float src[4], float dest[4])
{
    dest[0] = src[0];
    dest[1] = src[1];
    dest[2] = src[2];
    dest[3] = src[3];
}

float vec4fDot(float v0[4], float v1[4])
{
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] / v0[3] * v1[3];
}

void vec4fCross(float v0[4], float v1[4], float result[4])
{
    result[0] = v0[1] / v0[3] * v1[2] / v1[3] - v0[2] / v0[3] * v1[1] / v1[3];
    result[1] = v0[2] / v0[3] * v1[0] / v1[3] - v0[0] / v0[3] * v1[2] / v1[3];
    result[2] = v0[0] / v0[3] * v1[1] / v1[3] - v0[1] / v0[3] * v1[0] / v1[3];
    result[3] = 1.0f;  // ???
}

void vec4fBlend(float v0[4], float w0, float v1[4], float w1, float result[4])
{
    result[0] = v0[0] * w0 + v1[0] * w1;
    result[1] = v0[1] * w0 + v1[1] * w1;
    result[2] = v0[2] * w0 + v1[2] * w1;
    result[3] = 1.0f; // ??? 
}

void vec4fScale(float v[4], float w, float r[4])
{
    r[0] = v[0] * w;
    r[1] = v[1] * w;
    r[2] = v[2] * w;
    r[3] = v[3];
}

void vec4fNormalize(float v[4], float r[4])
{
    vec4fScale(v, 1.0 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]), r);
}

void vec4fReflect(float i[4], float n[4], float r[4])
{
    float dot;

    dot = vec4fDot(i, n);

    r[0] = i[0] - 2.0f * n[0] * dot;
    r[1] = i[1] - 2.0f * n[1] * dot;
    r[2] = i[2] - 2.0f * n[2] * dot;
    r[3] = i[3] - 2.0f * n[3] * dot; // This has to be wrong.
}

void mat4fMultVec4f(float m[16], float in[4], float out[4])
{
    int i;
    float t[4];

    for(i = 0; i < 4; i++)
	t[i] =
	    m[0 + i] * in[0] + 
	    m[4 + i] * in[1] + 
	    m[8 + i] * in[2] + 
	    m[12 + i] * in[3];
    memcpy(out, t, sizeof(out[0]) * 4);
}

void mat4fMultVec3fPt(float m[16], float in[3], float out[3])
{
    int i;
    float t[4];

    for(i = 0; i < 4; i++)
	t[i] =
	    m[0 + i] * in[0] + 
	    m[4 + i] * in[1] + 
	    m[8 + i] * in[2] + 
	    m[12 + i] * 1.0;

    t[0] /= t[3];
    t[1] /= t[3];
    t[2] /= t[3];

    memcpy(out, t, sizeof(out[0]) * 3);
}

void mat4fMultVec3fNm(float m[16], float in[3], float out[3])
{
    // multiply by inverse transpose
}

float mat4fDeterminant(float mat[16])
{
    return (mat[0] * mat[5] - mat[1] * mat[4]) *
        (mat[10] * mat[15] - mat[11] * mat[14]) + 
        (mat[2] * mat[4] - mat[0] * mat[6]) *
	(mat[9] * mat[15] - mat[11] * mat[13]) + 
        (mat[0] * mat[7] - mat[3] * mat[4]) *
	(mat[9] * mat[14] - mat[10] * mat[13]) + 
        (mat[1] * mat[6] - mat[2] * mat[5]) *
	(mat[8] * mat[15] - mat[11] * mat[12]) + 
        (mat[3] * mat[5] - mat[1] * mat[7]) *
	(mat[8] * mat[14] - mat[10] * mat[12]) + 
        (mat[2] * mat[7] - mat[3] * mat[6]) *
	(mat[8] * mat[13] - mat[9] * mat[12]);
}

float mat4fIdentity[16] =
{
    1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0,
};

#define EPSILON .00001

int mat4fInvert(float mat[16], float inv[16])
{
    int		i, rswap;
    float	det, div, swap;
    float	hold[16];

    memcpy(hold, mat, sizeof(mat[0]) * 16);
    memcpy(inv, mat4fIdentity, sizeof(mat[0]) * 16);
    det = mat4fDeterminant(mat);
    if(fabs(det) < EPSILON) /* singular? */
	return -1;

    rswap = 0;
    /* this loop isn't entered unless [0 + 0] > EPSILON and det > EPSILON,
	 so rswap wouldn't be 0, but I initialize so as not to get warned */
    if(fabs(hold[0]) < EPSILON)
    {
        if(fabs(hold[1]) > EPSILON)
            rswap = 1;
        else if(fabs(hold[2]) > EPSILON)
	    rswap = 2;
        else if(fabs(hold[3]) > EPSILON)
	    rswap = 3;

        for(i = 0; i < 4; i++)
	{
            swap = hold[i * 4 + 0];
            hold[i * 4 + 0] = hold[i * 4 + rswap];
            hold[i * 4 + rswap] = swap;

            swap = inv[i * 4 + 0];
            inv[i * 4 + 0] = inv[i * 4 + rswap];
            inv[i * 4 + rswap] = swap;
        }
    }
        
    div = hold[0];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 0] /= div;
        inv[i * 4 + 0] /= div;
    }

    div = hold[1];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 1] -= div * hold[i * 4 + 0];
        inv[i * 4 + 1] -= div * inv[i * 4 + 0];
    }
    div = hold[2];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 2] -= div * hold[i * 4 + 0];
        inv[i * 4 + 2] -= div * inv[i * 4 + 0];
    }
    div = hold[3];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 3] -= div * hold[i * 4 + 0];
        inv[i * 4 + 3] -= div * inv[i * 4 + 0];
    }

    if(fabs(hold[5]) < EPSILON){
        if(fabs(hold[6]) > EPSILON)
	    rswap = 2;
        else if(fabs(hold[7]) > EPSILON)
	    rswap = 3;

        for(i = 0; i < 4; i++)
	{
            swap = hold[i * 4 + 1];
            hold[i * 4 + 1] = hold[i * 4 + rswap];
            hold[i * 4 + rswap] = swap;

            swap = inv[i * 4 + 1];
            inv[i * 4 + 1] = inv[i * 4 + rswap];
            inv[i * 4 + rswap] = swap;
        }
    }

    div = hold[5];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 1] /= div;
        inv[i * 4 + 1] /= div;
    }

    div = hold[4];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 0] -= div * hold[i * 4 + 1];
        inv[i * 4 + 0] -= div * inv[i * 4 + 1];
    }
    div = hold[6];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 2] -= div * hold[i * 4 + 1];
        inv[i * 4 + 2] -= div * inv[i * 4 + 1];
    }
    div = hold[7];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 3] -= div * hold[i * 4 + 1];
        inv[i * 4 + 3] -= div * inv[i * 4 + 1];
    }

    if(fabs(hold[10]) < EPSILON){
        for(i = 0; i < 4; i++)
	{
            swap = hold[i * 4 + 2];
            hold[i * 4 + 2] = hold[i * 4 + 3];
            hold[i * 4 + 3] = swap;

            swap = inv[i * 4 + 2];
            inv[i * 4 + 2] = inv[i * 4 + 3];
            inv[i * 4 + 3] = swap;
        }
    }

    div = hold[10];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 2] /= div;
        inv[i * 4 + 2] /= div;
    }

    div = hold[8];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 0] -= div * hold[i * 4 + 2];
        inv[i * 4 + 0] -= div * inv[i * 4 + 2];
    }
    div = hold[9];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 1] -= div * hold[i * 4 + 2];
        inv[i * 4 + 1] -= div * inv[i * 4 + 2];
    }
    div = hold[11];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 3] -= div * hold[i * 4 + 2];
        inv[i * 4 + 3] -= div * inv[i * 4 + 2];
    }

    div = hold[15];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 3] /= div;
        inv[i * 4 + 3] /= div;
    }

    div = hold[12];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 0] -= div * hold[i * 4 + 3];
        inv[i * 4 + 0] -= div * inv[i * 4 + 3];
    }
    div = hold[13];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 1] -= div * hold[i * 4 + 3];
        inv[i * 4 + 1] -= div * inv[i * 4 + 3];
    }
    div = hold[14];
    for(i = 0; i < 4; i++)
    {
        hold[i * 4 + 2] -= div * hold[i * 4 + 3];
        inv[i * 4 + 2] -= div * inv[i * 4 + 3];
    }
    
    return 0;
}

static void mat4fMakeIdentity(float matrix[16])
{
    int i;
    for(i = 0; i < 16; i++)
	matrix[i] = 0.0f;
    matrix[0] = 1.0f;
    matrix[5] = 1.0f;
    matrix[10] = 1.0f;
    matrix[15] = 1.0f;
}

static void mat4fMakeTranslation(float x, float y, float z, float matrix[16])
{
    mat4fMakeIdentity(matrix);
    matrix[12] = x;
    matrix[13] = y;
    matrix[14] = z;
}

static void mat4fMakeScale(float x, float y, float z, float matrix[16])
{
    mat4fMakeIdentity(matrix);
    matrix[0] = x;
    matrix[5] = y;
    matrix[10] = z;
}

static void mat4fMakeRotation(float a, float x, float y, float z, float matrix[16])
{
    float c, s, t;

    c = (float)cos(a);
    s = (float)sin(a);
    t = 1.0f - c;

    matrix[0] = t * x * x + c;
    matrix[1] = t * x * y + s * z;
    matrix[2] = t * x * z - s * y;
    matrix[3] = 0;

    matrix[4] = t * x * y - s * z;
    matrix[5] = t * y * y + c;
    matrix[6] = t * y * z + s * x;
    matrix[7] = 0;

    matrix[8] = t * x * z + s * y;
    matrix[9] = t * y * z - s * x;
    matrix[10] = t * z * z + c;
    matrix[11] = 0;

    matrix[12] = 0;
    matrix[13] = 0;
    matrix[14] = 0;
    matrix[15] = 1;
}

static void mat4fGetRotation(float matrix[16], float rotation[4])
{
    float cosine;
    float sine;
    float d;

    cosine = (matrix[0] + matrix[5] + matrix[10] - 1.0f) / 2.0f;

    /* grantham 20000418 - I know this fixes the mysterious */
    /* NAN matrices, but I have no idea what the above number is supposed */
    /* to do, so I don't know why this happens */
    if(cosine > 1.0){
#if defined(DEBUG)
	fprintf(stderr, "XXX acos of greater than 1! (clamped)\n");
#endif // DEBUG
	cosine = 1.0;
    }
    if(cosine < -1.0){
#if defined(DEBUG)
	fprintf(stderr, "XXX acos of less than -1! (clamped)\n");
#endif // DEBUG
	cosine = -1.0;
    }

    rotation[0] = (float)acos(cosine);

#if defined(DEBUG)
    if(rotation[0] != rotation[0]) /* isNAN */
	abort();
#endif // DEBUG

    sine = (float)sin(rotation[0]);
    rotation[1] = (matrix[6] - matrix[9]);
    rotation[2] = (matrix[8] - matrix[2]);
    rotation[3] = (matrix[1] - matrix[4]);
    d = sqrt(rotation[1] * rotation[1] + rotation[2] * rotation[2] +
	rotation[3] * rotation[3]);
    rotation[1] /= d;
    rotation[2] /= d;
    rotation[3] /= d;
}

static void mat4fMultMat4f(float m1[16], float m2[16], float r[16])
{
    float t[16];
    int i, j;

    for(j = 0; j < 4; j++)
	for(i = 0; i < 4; i++)
           t[i * 4 + j] = m1[i * 4 + 0] * m2[0 * 4 + j] +
	       m1[i * 4 + 1] * m2[1 * 4 + j] +
	       m1[i * 4 + 2] * m2[2 * 4 + j] +
	       m1[i * 4 + 3] * m2[3 * 4 + j];

    memcpy(r, t, sizeof(t));
}

static void rot4MultRot4(float rotation1[4], float rotation2[4], float result[4])
{
    float matrix1[16];
    float matrix2[16];
    float matrix3[16];
    float dist;

    mat4fMakeRotation(rotation1[0], rotation1[1], rotation1[2], rotation1[3],
        matrix1);
    mat4fMakeRotation(rotation2[0], rotation2[1], rotation2[2], rotation2[3],
        matrix2);
    mat4fMultMat4f(matrix1, matrix2, matrix3);
    mat4fGetRotation(matrix3, result);

    dist = (float)sqrt(result[1] * result[1] + result[2] * result[2] +
        result[3] * result[3]);

#if defined(DEBUG)
    if(result[0] != result[0]) /* isNAN */
	abort();
#endif // DEBUG

    result[1] /= dist;
    result[2] /= dist;
    result[3] /= dist;
}


#define MAXX 0
#define MINX 1
#define MAXY 2
#define MINY 3
#define MAXZ 4
#define MINZ 5

void boxInit(float b[6])
{
    b[MAXX] = b[MAXY] = b[MAXZ] = -FLT_MAX;
    b[MINX] = b[MINY] = b[MINZ] = FLT_MAX;
}

void boxExtendByPoint(float b[6], float x, float y, float z)
{
    if(x < b[MINX]) b[MINX] = x;
    if(x > b[MAXX]) b[MAXX] = x;
    if(y < b[MINY]) b[MINY] = y;
    if(y > b[MAXY]) b[MAXY] = y;
    if(z < b[MINZ]) b[MINZ] = z;
    if(z > b[MAXZ]) b[MAXZ] = z;
}

void boxExtendBySphere(float b[6], float x, float y, float z, float r)
{
    if(x - r < b[MINX]) b[MINX] = x - r;
    if(x + r > b[MAXX]) b[MAXX] = x + r;
    if(y - r < b[MINY]) b[MINY] = y - r;
    if(y + r > b[MAXY]) b[MAXY] = y + r;
    if(z - r < b[MINZ]) b[MINZ] = z - r;
    if(z + r > b[MAXZ]) b[MAXZ] = z + r;
}


#include <math.h>

enum xformMode {
    XFORM_MODE_TRACKBALL,	/* rotate in direction of mouse motion */
    XFORM_MODE_ROLL,		/* rotate around Z axis */
    XFORM_MODE_SCROLL,		/* translate in X-Y */
    XFORM_MODE_DOLLY		/* translate in Z */
};

typedef struct { 
    enum xformMode mode;
    float matrix[16];

    float frame[16];		/* The coordinate frame for this transform */

    float worldX[3];		/* world X axis in this coordinate space */
    float worldY[3];		/* world Y axis in this coordinate space */	
    float worldZ[3];		/* world Z axis in this coordinate space */

    float referenceSize;	/* used to calculate translations */
    float motionScale;		/* for dynamic scaling, etc */

    float rotation[4];		/* radians, x, y, z, like OpenGL */
    float translation[3];
    float scale[3];		/* scaled around center. */
    float center[3];		/* ignore by setting to <0,0,0> */
} Transform;

/* Win32 math.h doesn't define M_PI. */
#if defined(_WIN32)
#if !defined(M_PI)
#define M_PI 3.14159265
#endif /* !defined(M_PI) */
#endif /* defined(WIN32) */

/* internal use */

static void dragToRotation(float dx, float dy, float *rotation)
{
    float dist;

    /* XXX grantham 990825 - this "dist" doesn't make me confident. */
    /* but I put in the *10000 to decrease chance of underflow  (???) */
    dist = sqrt(dx * 10000 * dx * 10000 + dy * 10000 * dy * 10000) / 10000;
    /* dist = sqrt(dx * dx + dy * dy); */
    rotation[0] = (float) M_PI * dist;
    rotation[1] = (float) dy / dist;
    rotation[2] = (float) dx / dist;
    rotation[3] = 0.0f;
}

static void calcViewMatrix(float viewRotation[4], float viewOffset[3],
    float objCenter[3], float objScale[3], float viewMatrix[16])
{
    float tmp[16];

    /* These could be generated with OpenGL matrix functions */
    mat4fMakeIdentity(viewMatrix);
    mat4fMakeTranslation(viewOffset[0], viewOffset[1], viewOffset[2], tmp);
    mat4fMultMat4f(tmp, viewMatrix, viewMatrix);
    mat4fMakeRotation(viewRotation[0], viewRotation[1], viewRotation[2],
        viewRotation[3], tmp);
    mat4fMultMat4f(tmp, viewMatrix, viewMatrix);
    mat4fMakeScale(objScale[0], objScale[1], objScale[2], tmp);
    mat4fMultMat4f(tmp, viewMatrix, viewMatrix);
    mat4fMakeTranslation(-objCenter[0], -objCenter[1], -objCenter[2], tmp);
    mat4fMultMat4f(tmp, viewMatrix, viewMatrix);
}

/* external API */

void xformCalcMatrix(Transform *xform)
{
    calcViewMatrix(xform->rotation, xform->translation, xform->center,
        xform->scale, xform->matrix);
}

void xformSetFrame(Transform *xform, float *frame)
{
    float i[16];
    float origin[3];

    memcpy(xform->frame, frame, sizeof(float) * 16); 

    mat4fInvert(xform->frame, i);
    mat4fMultVec3fPt(i, vec3fOrigin, origin);
    mat4fMultVec3fPt(i, vec3fXAxis, xform->worldX);
    mat4fMultVec3fPt(i, vec3fYAxis, xform->worldY);
    mat4fMultVec3fPt(i, vec3fZAxis, xform->worldZ);
    vec3fBlend(xform->worldX, 1.0f, origin, -1.0f, xform->worldX);
    vec3fBlend(xform->worldY, 1.0f, origin, -1.0f, xform->worldY);
    vec3fBlend(xform->worldZ, 1.0f, origin, -1.0f, xform->worldZ);
    /* figure out new X, Y, Z */
}

void xformMotion(Transform *xform, float dx, float dy)
{
    switch(xform->mode) {

	case XFORM_MODE_TRACKBALL:
	    if(dx != 0 || dy != 0) {
		float localRotation[4];
		float worldRotation[4];
		dragToRotation(dx, dy, worldRotation);
		localRotation[0] = worldRotation[0];
		vec3fScale(xform->worldX, worldRotation[1], localRotation + 1);
		vec3fBlend(localRotation + 1, 1.0f, xform->worldY, worldRotation[2], localRotation + 1);
		rot4MultRot4(xform->rotation, localRotation, xform->rotation);
	    }
	    break;

	case XFORM_MODE_ROLL:
	    {
		float rotation[4];
		rotation[0] = M_PI * 2 * -dy;
		rotation[1] = 0.0f;
		rotation[2] = 0.0f;
		rotation[3] = 1.0f;
		rot4MultRot4(xform->rotation, rotation, xform->rotation);
	    }
	    break;

	case XFORM_MODE_SCROLL:
	    vec3fBlend(xform->translation, 1.0f, xform->worldX, dx * xform->referenceSize * xform->motionScale, xform->translation);
	    vec3fBlend(xform->translation, 1.0f, xform->worldY, -dy * xform->referenceSize * xform->motionScale, xform->translation);
	    break;

	case XFORM_MODE_DOLLY:
	    vec3fBlend(xform->translation, 1.0f, xform->worldZ, dy * xform->referenceSize * xform->motionScale, xform->translation);
	    break;
    }

    calcViewMatrix(xform->rotation, xform->translation, xform->center,
        xform->scale, xform->matrix);
}

void xformInitialize(Transform *xform)
{
    vec3fCopy(vec3fXAxis, xform->worldX);
    vec3fCopy(vec3fYAxis, xform->worldY);
    vec3fCopy(vec3fZAxis, xform->worldZ);

    xform->center[0] = 0.0f;
    xform->center[1] = 0.0f;
    xform->center[2] = 0.0f;

    xform->scale[0] = 1.0f;
    xform->scale[1] = 1.0f;
    xform->scale[2] = 1.0f;

    /* diagonal of 2-high cube */
    xform->referenceSize = 3.465;

    xform->motionScale = 1.0;

    xform->translation[0] = 0.0f;
    xform->translation[1] = 0.0f;
    xform->translation[2] = 0.0f;

    xform->rotation[0] = 0.0;
    xform->rotation[1] = 1.0;
    xform->rotation[2] = 0.0;
    xform->rotation[3] = 0.0;

    calcViewMatrix(xform->rotation, xform->translation, xform->center,
        xform->scale, xform->matrix);

    xform->mode = XFORM_MODE_TRACKBALL;
}


void xformInitializeViewFromBox(Transform *xform, float box[6], float fov)
{
    vec3fCopy(vec3fXAxis, xform->worldX);
    vec3fCopy(vec3fYAxis, xform->worldY);
    vec3fCopy(vec3fZAxis, xform->worldZ);

    xform->center[0] = (box[MINX] + box[MAXX]) / 2.0f;
    xform->center[1] = (box[MINY] + box[MAXY]) / 2.0f;
    xform->center[2] = (box[MINZ] + box[MAXZ]) / 2.0f;

    xform->scale[0] = 1.0f;
    xform->scale[1] = 1.0f;
    xform->scale[2] = 1.0f;

    if(box[MAXX] - box[MINX] > box[MAXY] - box[MINY] &&
	box[MAXX] - box[MINX] > box[MAXZ] - box[MINZ])
        xform->referenceSize = box[MAXX] - box[MINX];
    else if(box[MAXY] - box[MINY] > box[MAXZ] - box[MINZ])
        xform->referenceSize = box[MAXY] - box[MINY];
    else
        xform->referenceSize = box[MAXZ] - box[MINZ];

    xform->motionScale = 1.0;

    xform->translation[0] = 0.0f;
    xform->translation[1] = 0.0f;
    xform->translation[2] = -xform->referenceSize / cos(fov / 2.0);

    xform->rotation[0] = 0.0;
    xform->rotation[1] = 1.0;
    xform->rotation[2] = 0.0;
    xform->rotation[3] = 0.0;

    calcViewMatrix(xform->rotation, xform->translation, xform->center,
        xform->scale, xform->matrix);

    xform->mode = XFORM_MODE_TRACKBALL;
}

struct String {
    float x, y;
    float r, g, b;
    string s; 
};

struct Line {
    float x1, y1;
    float x2, y2;
    float r, g, b;
    float width;
};

struct Point {
    float x, y;
    float r, g, b;
    float size;
};

vector<String> strings;
vector<Line> lines;
vector<Point> points;
/* icons later */

void readArt (FILE *fp)
{
    char input[512];

    while(fgets(input, sizeof(input) - 1, fp) != NULL) {

	input[strlen(input) - 1] = '\0';

	if(strncmp("string ", input, 7) == 0) {

            int soFar;
            float x, y, r, g, b;
            sscanf(input + 7, "%f %f %f %f %f %n", &x, &y,
                &r, &g, &b, &soFar);
            string s(input + 7 + soFar);
            String str = {x, y, r, g, b, s};
            strings.push_back(str);

	} else if(strncmp("line ", input, 5) == 0) {

            float x1, x2, y1, y2, r, g, b, width;
            sscanf(input + 5, "%f %f %f %f %f %f %f %f",
                &x1, &y1, &x2, &y2, &r, &g, &b, &width);
            Line line = {x1, y1, x2, y2, r, g, b, width};
            lines.push_back(line);

	} else if(strncmp("point ", input, 6) == 0) {

            float x, y, r, g, b, size;
            sscanf(input + 5, "%f %f %f %f %f %f",
                &x, &y, &r, &g, &b, &size);
            Point point = {x, y, r, g, b, size};
            points.push_back(point);
	}
    }
}

float sceneBox[6];

void makeBoxes(void)
{
    boxInit(sceneBox);

    for(auto it = strings.begin(); it != strings.end(); it++) {
	String& s = *it;
	boxExtendByPoint(sceneBox, s.x, s.y, 0.0);
    }

    for(auto it = lines.begin(); it != lines.end(); it++) {
	Line& l = *it;
	boxExtendByPoint(sceneBox, l.x1, l.y1, 0.0);
	boxExtendByPoint(sceneBox, l.x2, l.y2, 0.0);
    }

    for(auto it = points.begin(); it != points.end(); it++) {
	Point& p = *it;
	boxExtendByPoint(sceneBox, p.x, p.y, 0.0);
    }
}

void
drawString(void *font, char *str)
{
    int len, i;

    len = (int) strlen(str);
    for (i = 0; i < len; i++) {
        glutBitmapCharacter(font, str[i]);
    }
}

int                  winWidth, winHeight;

void glstring(void *font, float x, float y, const char *fmt, ...)
{
    va_list args;
    char dummy[512];

    glRasterPos2f(x, y);

    va_start(args, fmt);
    vsprintf(dummy, fmt, args);
    va_end(args);

    drawString(font, dummy);
}

void glprintf(void *font, int x, int y, char *fmt, ...)
{
    va_list args;
    char dummy[512];

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    glOrtho(0, winWidth, 0, winHeight, -1, 1);
    glRasterPos2f(x, y);

    va_start(args, fmt);
    vsprintf(dummy, fmt, args);
    va_end(args);

    drawString(font, dummy);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
}

void drawStrings(void)
{
    for(auto it = strings.begin(); it != strings.end(); it++) {
	String& s = *it;

	glColor3f(s.r, s.g, s.b);
	glstring(GLUT_BITMAP_9_BY_15, s.x, s.y, " %s", s.s.c_str());
    }
}

void drawLines(void)
{
    for(auto it = lines.begin(); it != lines.end(); it++) {
	Line& l = *it;

	glColor3f(l.r, l.g, l.b);
	glLineWidth(l.width);
	glBegin(GL_LINES);
	glVertex2f(l.x1, l.y1);
	glVertex2f(l.x2, l.y2);
	glEnd();
    }
}

void drawPoints(void)
{
    for(auto it = points.begin(); it != points.end(); it++) {
	Point& p = *it;

	glColor3f(p.r, p.g, p.b);
	glPointSize(p.size);
	glBegin(GL_POINTS);
	glVertex2f(p.x, p.y);
	glEnd();
    }
}

Transform            mainTransform;

float                aspectRatio = 1;

void init(void)
{
    xformInitializeViewFromBox(&mainTransform, sceneBox, .57595f);

    glClearColor(.2, .2, .4, 1);
}

void redraw(void)
{ 
    float nearClip, farClip;

    nearClip = - mainTransform.translation[2] - mainTransform.referenceSize;
    farClip = - mainTransform.translation[2] + mainTransform.referenceSize;
    if(nearClip < 0.1 * mainTransform.referenceSize)
        nearClip = 0.1 * mainTransform.referenceSize;
    if(farClip < 0.2 * mainTransform.referenceSize)
        nearClip = 0.2 * mainTransform.referenceSize;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-.66 / aspectRatio * nearClip, .66 / aspectRatio * nearClip,
        -.66 * nearClip, .66 * nearClip,
        nearClip, farClip);

    glMatrixMode(GL_MODELVIEW);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();

    glMultMatrixf((float *)mainTransform.matrix);

    drawStrings();
    drawLines();
    drawPoints();

    glPopMatrix();

    glutSwapBuffers();
}

void keyup(unsigned char key, int x, int y)
{
    switch(key) {
    case 'f':
        break;
    }
}

void keyboard(unsigned char key, int x, int y)
{
    switch(key)
    {
        case 'r':
            mainTransform.mode = XFORM_MODE_TRACKBALL;
            break;

        case 'o':
            mainTransform.mode = XFORM_MODE_ROLL;
            break;

        case 'x':
            mainTransform.mode = XFORM_MODE_SCROLL;
            break;

        case 'z':
            mainTransform.mode = XFORM_MODE_DOLLY;
            break;

        case 'q': case 'Q': case '\033':
            exit(0);
            break;

        case 's':
            {
                static unsigned char *pixels;
                static int pixelsSize;
                int viewport[4];
                int i;
                FILE *fp;

                if((fp = fopen("screenshot.ppm", "wb")) == NULL) {
                    fprintf(stderr, "snapshot: couldn't open "
                        "\"screenshot.ppm\".\n");
                    goto noshot;
                }
                glGetIntegerv(GL_VIEWPORT, viewport);
                if(pixelsSize < viewport[2] * viewport[3] * 3) {
                    pixelsSize = viewport[2] * viewport[3] * 3;
                    pixels = (unsigned char *)realloc(pixels, pixelsSize);
                    if(pixels == NULL) {
                        fprintf(stderr, "snapshot: couldn't allocate %u "
                            "bytes for screenshot.\n", pixelsSize);
                        goto noshot;
                    }
                }
                glPushAttrib(GL_CLIENT_PIXEL_STORE_BIT);
                glPixelStorei(GL_PACK_ALIGNMENT, 1);
                glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3],
                    GL_RGB, GL_UNSIGNED_BYTE, pixels);
                glPopAttrib();
                fprintf(fp, "P6 %d %d 255\n", viewport[2], viewport[3]);
                for(i = viewport[3] - 1; i >= 0; i--)
                    fwrite(pixels + viewport[2] * 3 * i, 3, viewport[2], fp);
                fclose(fp);
        noshot:;
            }
            break;


        default:
            printf("I don't know what do do about key '%c'\n", key);
            break;
    }
}


static int ox, oy;


void button(int b, int state, int x, int y)
{
    if(state == 1) {
    } else {
        ox = x;
        oy = y;
    }
}


void motion(int x, int y)
{
    int dx, dy;

    dx = x - ox;
    dy = y - oy;

    ox = x;
    oy = y;

    xformMotion(&mainTransform, dx / 512.0, dy / 512.0);
}


void reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    winWidth = width;
    winHeight = height;
    aspectRatio = height / (float) width;
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    int argsused;
    FILE *in;
    char *progname = argv[0];

    glutInit(&argc, argv);
    glutInitWindowSize(512, 512);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutCreateWindow("simple primitive viewer (again)");
    glutDisplayFunc(redraw);
    glutIdleFunc(redraw);
    glutKeyboardUpFunc(keyup);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(motion);
    glutMouseFunc(button);
    glutReshapeFunc(reshape);

    argc--;
    argv++;
    while(argc > 0 && strcmp(argv[0], "-") != 0 && argv[0][0] == '-') {
        argsused = 1;
        if(strcmp(argv[0], "-h") == 0) {
            fprintf(stderr, "usage: %s [options] [--] filename\n", progname);
            fprintf(stderr, "use \"-\" for stdin.\n");
            fprintf(stderr, "\te.g. %s -\n", progname);
            fprintf(stderr, "options:\n");
        } else {
            fprintf(stderr, "didn't understand command line parameter "
                "\"%s\"\n", argv[0]);
            exit(1);
        }
        argc -= argsused;
        argv += argsused;
    }
    if(argc < 1) {
        fprintf(stderr, "Expected filename parameter.\n");
        exit(1);
    }

    if(strcmp(argv[0], "-") == 0) {
        in = stdin;
    } else {
        in = fopen(argv[0], "r");
        if(in == NULL) {
             fprintf(stderr, "Couldn't open \"%s\" for reading\n", argv[0]);
             exit(1);
        }
    }

    readArt(in);
    makeBoxes();


    init();
    glutMainLoop();
    return(0);
}

/* vi:ts=8 sw=4
 */
