///////////////////////////////////////////////////////////////////////////
//
// NAME
//  Convert.cpp -- convert between image types, copy images, select bands
//
// DESCRIPTION
//  See Convert.h for a full description
//
// SEE ALSO
//  Convert.h           full description
//
// Copyright © Richard Szeliski, 2001.
// See Copyright.h for more details
//
///////////////////////////////////////////////////////////////////////////

#include "Image.h"
#include "Error.h"
#include "Convert.h"

//
// Type conversion utilities
//

template <class T1, class T2>
extern void ScaleAndOffsetLine(T1* src, T2* dst, int n,
                        float scale, float offset,
                        T2 minVal, T2 maxVal)
{
    // This routine does NOT round values when converting from float to int
    const bool scaleOffset = (scale != 1.0f) || (offset != 0.0f);
    const bool clip = (minVal < maxVal);

    if (scaleOffset)
        for (int i = 0; i < n; i++)
        {
            float val = src[i] * scale + offset;
            if (clip)
                val = __min(__max(val, minVal), maxVal);
            dst[i] = (T2) val;
        }
    else if (clip)
        for (int i = 0; i < n; i++)
        {
            dst[i] = (T2) __min(__max(src[i], minVal), maxVal);
        }
    else if (typeid(T1) == typeid(T2))
        memcpy(dst, src, n*sizeof(T2));
    else
        for (int i = 0; i < n; i++)
        {
            dst[i] = (T2) src[i];
        }
}

template <class T1, class T2>
extern void ScaleAndOffset(CImageOf<T1>& src, CImageOf<T2>& dst, float scale, float offset)
{
    // Convert between images of same shape but diffent types,
    //  and optionally scale and/or offset the pixel values
    CShape sShape = src.Shape();
    CShape dShape = dst.Shape();

    // Make sure the shapes (ignoring bands) are compatible
    if (sShape != dShape)
        dst.ReAllocate(sShape);

    // Determine if clipping is required
    T2 minVal = dst.MinVal();
    T2 maxVal = dst.MaxVal();

    // if (minVal <= src.MinVal() && maxVal >= src.MaxVal())
    // changed DS 4/2003 - take into consideration scale and offset
    // removed DS 1/2004 - doesn't work for negative scales - safer to always clip
    //if (src.MinVal() * scale + offset >= minVal && 
    //    src.MaxVal() * scale + offset <= maxVal)
    //    minVal = maxVal = 0;

    // Process each row
    for (int y = 0; y < sShape.height; y++)
    {
        int n = sShape.width * sShape.nBands;
        ScaleAndOffsetLine(&src.Pixel(0, y, 0), &dst.Pixel(0, y, 0),
                           n, scale, offset, minVal, maxVal);
    }
}

template <class T>
extern CImageOf<T> ConvertToRGBA(CImageOf<T> src)
{
    // Check if already RGBA
    CShape sShape = src.Shape();
    if (sShape.nBands == 4 && src.alphaChannel == 3)
        return src;

    // Make sure the source is a gray image
    if (sShape.nBands != 1)
        throw CError("ConvertToRGBA: can only convert from 1-band (gray) image");

    // Allocate the new image
    CShape dShape(sShape.width, sShape.height, 4);
    CImageOf<T> dst(dShape);

    // Process each row
    int aC = dst.alphaChannel;
    for (int y = 0; y < sShape.height; y++)
    {
        T* srcP = &src.Pixel(0, y, 0);
        T* dstP = &dst.Pixel(0, y, 0);
        for (int x = 0; x < sShape.width; x++, srcP++)
            for (int b = 0; b < dShape.nBands; b++, dstP++)
                *dstP = (b == aC) ? 255 : *srcP;
    }
    return dst;
}

template <class T>
extern CImageOf<T> ConvertToGray(CImageOf<T> src)
{
    // Check if already gray
    CShape sShape = src.Shape();
    if (sShape.nBands == 1)
        return src;

    // Make sure the source is a color image
    if (sShape.nBands != 4 || src.alphaChannel != 3)
        throw CError("ConvertToGray: can only convert from 4-band (RGBA) image");

    // Allocate the new image
    CShape dShape(sShape.width, sShape.height, 1);
    CImageOf<T> dst(dShape);

    // Process each row
    T minVal = dst.MinVal();
    T maxVal = dst.MaxVal();
    for (int y = 0; y < sShape.height; y++)
    {
        T* srcP = &src.Pixel(0, y, 0);
        T* dstP = &dst.Pixel(0, y, 0);
        for (int x = 0; x < sShape.width; x++, srcP += 4, dstP++)
        {
            RGBA<T>& p = *(RGBA<T> *) srcP;
            // OLD FORMULA: float Y = (float)(0.212671 * p.R + 0.715160 * p.G + 0.072169 * p.B);
	    // Changed to commonly used formula 6/4/07 DS
            float Y = (float)(0.299 * p.R + 0.587 * p.G + 0.114 * p.B);
            *dstP = (T) __min(maxVal, __max(minVal, Y));
        }
    }
    return dst;
}

template <class T>
extern void BandSelect(CImageOf<T>& src, CImageOf<T>& dst, int sBand, int dBand)
{
    // Convert between images of same type but different # of bands
    CShape sShape = src.Shape();
    CShape dShape = dst.Shape();

    // Make sure the shapes (ignoring bands) are compatible
    if (! sShape.SameIgnoringNBands(dShape) || dShape.nBands == 0)
    {
        dShape.width  = sShape.width;
        dShape.height = sShape.height;
        dShape.nBands = (dShape.nBands) ? dShape.nBands : 1;
        dst.ReAllocate(dShape);
    }

    // Check the bands are valid
    int sB = sShape.nBands;
    int dB = dShape.nBands;
    if (sBand < 0 || sBand >= sB)
        throw CError("BandSelect: source band %d is invalid", sBand);
    if (dBand < 0 || dBand >= dB)
        throw CError("BandSelect: destination band %d is invalid", dBand);

    // Process each row
    for (int y = 0; y < sShape.height; y++)
    {
        T* srcP = &src.Pixel(0, y, 0);
        T* dstP = &dst.Pixel(0, y, 0);
        for (int x = 0; x < sShape.width; x++, srcP += sB, dstP += dB)
            dstP[dBand] = srcP[sBand];
    }
}

//
// Force instantiation for the types we care about (uchar, int, float)
//

// for some reason the following function doesn't work under gcc 4 with -O2:
/*
template <class T1>
void CopyPixelsInstantiate(CImageOf<T1>& s1)
{
    CByteImage  b2;
    CIntImage   i2;
    CFloatImage f2;
    CopyPixels(s1, b2);
    CopyPixels(s1, i2);
    CopyPixels(s1, f2);
}
*/
// ... so do it like this instead:
template void ScaleAndOffset(CByteImage&  src, CByteImage&  dst, float s, float o);
template void ScaleAndOffset(CByteImage&  src, CIntImage&   dst, float s, float o);
template void ScaleAndOffset(CByteImage&  src, CFloatImage& dst, float s, float o);
template void ScaleAndOffset(CIntImage&   src, CByteImage&  dst, float s, float o);
template void ScaleAndOffset(CIntImage&   src, CIntImage&   dst, float s, float o);
template void ScaleAndOffset(CIntImage&   src, CFloatImage& dst, float s, float o);
template void ScaleAndOffset(CFloatImage& src, CByteImage&  dst, float s, float o);
template void ScaleAndOffset(CFloatImage& src, CIntImage&   dst, float s, float o);
template void ScaleAndOffset(CFloatImage& src, CFloatImage& dst, float s, float o);

// also need (for Convolve) at least this:
template void ScaleAndOffsetLine(float* src, float* dst, int n, float scale, float offset, float minVal, float maxVal);

// same here:

template CByteImage ConvertToGray(CByteImage src);
template CIntImage ConvertToGray(CIntImage src);
template CFloatImage ConvertToGray(CFloatImage src);

template CByteImage ConvertToRGBA(CByteImage src);
template CIntImage ConvertToRGBA(CIntImage src);
template CFloatImage ConvertToRGBA(CFloatImage src);

template void BandSelect(CByteImage&  src, CByteImage&  dst, int sBand, int dBand);
template void BandSelect(CIntImage&   src, CIntImage&   dst, int sBand, int dBand);
template void BandSelect(CFloatImage& src, CFloatImage& dst, int sBand, int dBand);

/*
template <class T>
int InstantiateConvert(CImageOf<T> src)
{
    CopyPixelsInstantiate(src);
    CImageOf<T> r1 = ConvertToRGBA(src);
    CImageOf<T> r2 = ConvertToGray(src);
    BandSelect(r1, r2, 0, 0);
    return (int)r1.Pixel(0, 0, 0);
}

void InstantiateAllConverts(void)
{
    InstantiateConvert(CByteImage());
    InstantiateConvert(CIntImage());
    InstantiateConvert(CFloatImage());
}
*/
