/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkScalarsToColors.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkScalarsToColors.h"
#include "vtkTemplateAliasMacro.h"
#include "vtkUnsignedCharArray.h"
#include "vtkObjectFactory.h"

#include <math.h>

vtkStandardNewMacro(vtkScalarsToColors);

//----------------------------------------------------------------------------
vtkScalarsToColors::vtkScalarsToColors()
{
  this->Alpha = 1.0;
  this->VectorComponent = 0;
  this->VectorMode = vtkScalarsToColors::COMPONENT;
  this->UseMagnitude = 0;

  // only used in this class, not used in subclasses
  this->InputRange[0] = 0.0;
  this->InputRange[1] = 255.0;
}

//----------------------------------------------------------------------------
// Description:
// Return true if all of the values defining the mapping have an opacity
// equal to 1. Default implementation return true.
int vtkScalarsToColors::IsOpaque()
{
  return 1;
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::SetVectorModeToComponent()
{
  this->SetVectorMode(vtkScalarsToColors::COMPONENT);
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::SetVectorModeToMagnitude()
{
  this->SetVectorMode(vtkScalarsToColors::MAGNITUDE);
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::SetVectorModeToRGBColors()
{
  this->SetVectorMode(vtkScalarsToColors::RGBCOLORS);
}

//----------------------------------------------------------------------------
// do not use SetMacro() because we do not want the table to rebuild.
void vtkScalarsToColors::SetAlpha(double alpha)
{
  this->Alpha = (alpha < 0.0 ? 0.0 : (alpha > 1.0 ? 1.0 : alpha));
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::SetRange(double minval, double maxval)
{
  if (this->InputRange[0] != minval ||
      this->InputRange[1] != maxval)
    {
    this->InputRange[0] = minval;
    this->InputRange[1] = maxval;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
double *vtkScalarsToColors::GetRange()
{
  return this->InputRange;
}

//----------------------------------------------------------------------------
vtkIdType vtkScalarsToColors::GetNumberOfAvailableColors()
{
  // return total possible RGB colors
  return 256*256*256;
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::DeepCopy(vtkScalarsToColors *obj)
{
  if (obj)
    {
    this->Alpha = obj->Alpha;
    this->UseMagnitude = obj->UseMagnitude;
    this->VectorMode = obj->VectorMode;
    this->VectorComponent = obj->VectorComponent;
    this->InputRange[0] = obj->InputRange[0];
    this->InputRange[1] = obj->InputRange[1];
    }
}

//----------------------------------------------------------------------------
inline void vtkScalarsToColorsComputeShiftScale(
  vtkScalarsToColors *self, double &shift, double &scale)
{
  static double minscale = -1e17;
  static double maxscale = 1e17;

  double *range = self->GetRange();
  shift = -range[0];
  scale = range[1] - range[0];
  if (scale*scale > 1e-30)
    {
    scale = 1.0/scale;
    }
  else
    {
    scale = (scale < 0 ? minscale : maxscale);
    }
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::GetColor(double v, double rgb[3])
{
  static double minval = 0.0;
  static double maxval = 1.0;

  double shift, scale;
  vtkScalarsToColorsComputeShiftScale(this, shift, scale);

  double val = (v + shift)*scale;
  val = (val > minval ? val : minval);
  val = (val < maxval ? val : maxval);

  rgb[0] = val;
  rgb[1] = val;
  rgb[2] = val;
}

//----------------------------------------------------------------------------
double vtkScalarsToColors::GetOpacity(double vtkNotUsed(v))
{
  return 1.0;
}

//----------------------------------------------------------------------------
unsigned char *vtkScalarsToColors::MapValue(double v)
{
  double rgb[3];

  this->GetColor(v, rgb);
  double alpha = this->GetOpacity(v);

  this->RGBABytes[0] = static_cast<unsigned char>(rgb[0]*255 + 0.5);
  this->RGBABytes[1] = static_cast<unsigned char>(rgb[1]*255 + 0.5);
  this->RGBABytes[2] = static_cast<unsigned char>(rgb[2]*255 + 0.5);
  this->RGBABytes[3] = static_cast<unsigned char>(alpha*255 + 0.5);

  return this->RGBABytes;
}

//----------------------------------------------------------------------------
vtkUnsignedCharArray *vtkScalarsToColors::MapScalars(vtkDataArray *scalars,
                                                     int colorMode, int comp)
{
  int numberOfComponents = scalars->GetNumberOfComponents();
  vtkUnsignedCharArray *newColors;
  vtkUnsignedCharArray *colors;

  // map scalars through lookup table only if needed
  if ( colorMode == VTK_COLOR_MODE_DEFAULT && 
       (colors=vtkUnsignedCharArray::SafeDownCast(scalars)) != NULL )
    {
    newColors = this->
      ConvertUnsignedCharToRGBA(colors, colors->GetNumberOfComponents(),
                                scalars->GetNumberOfTuples());
    }
  else
    {
    newColors = vtkUnsignedCharArray::New();
    newColors->SetNumberOfComponents(4);
    newColors->SetNumberOfTuples(scalars->GetNumberOfTuples());

    // If mapper did not specify a component, use the VectorMode
    if (comp < 0 && numberOfComponents > 1)
      {
      this->MapVectorsThroughTable(scalars->GetVoidPointer(0),
                                   newColors->GetPointer(0),
                                   scalars->GetDataType(),
                                   scalars->GetNumberOfTuples(),
                                   scalars->GetNumberOfComponents(),
                                   VTK_RGBA);
      }
    else
      {
      if (comp < 0)
        {
        comp = 0;
        }
      if (comp >= numberOfComponents)
        {
        comp = numberOfComponents - 1;
        }

      // Map the scalars to colors
      this->MapScalarsThroughTable(scalars->GetVoidPointer(comp),
                                   newColors->GetPointer(0),
                                   scalars->GetDataType(),
                                   scalars->GetNumberOfTuples(),
                                   scalars->GetNumberOfComponents(),
                                   VTK_RGBA);
      }
    }

  return newColors;
}

//----------------------------------------------------------------------------
// Map a set of vector values through the table
void vtkScalarsToColors::MapVectorsThroughTable(
  void *input, unsigned char *output, int scalarType,
  int numValues, int inComponents, int outputFormat)
{
  if (outputFormat < VTK_LUMINANCE || outputFormat > VTK_RGBA)
    {
    vtkErrorMacro(<< "MapVectorsThroughTable: unrecognized color format");
    return;
    }

  this->UseMagnitude = 0;

  switch(this->GetVectorMode())
    {
    case vtkScalarsToColors::COMPONENT:
      {
      int scalarSize = vtkDataArray::GetDataTypeSize(scalarType);
      int vectorComponent = this->GetVectorComponent();
      if (vectorComponent < 0)
        {
        vectorComponent = 0;
        }
      if (vectorComponent >= inComponents)
        {
        vectorComponent = inComponents - 1;
        }

      this->MapScalarsThroughTable2(
        static_cast<unsigned char *>(input) + vectorComponent*scalarSize,
        output, scalarType, numValues, inComponents, outputFormat);
      }
      break;

    case vtkScalarsToColors::MAGNITUDE:
      {
      this->UseMagnitude = 1;
      this->MapScalarsThroughTable2(
        input, output, scalarType, numValues, inComponents, outputFormat);
      }
      break;

    case vtkScalarsToColors::RGBCOLORS:
      {
      this->MapColorsToColors(
        input, output, scalarType, numValues, inComponents, outputFormat);
      }
      break;
   }
}

//----------------------------------------------------------------------------
// Map a set of scalar values through the table
void vtkScalarsToColors::MapScalarsThroughTable(vtkDataArray *scalars, 
                                                unsigned char *output,
                                                int outputFormat)
{
  if (outputFormat < VTK_LUMINANCE || outputFormat > VTK_RGBA)
    {
    vtkErrorMacro(<< "MapScalarsThroughTable: unrecognized color format");
    return;
    }

  this->MapScalarsThroughTable(scalars->GetVoidPointer(0),
                               output,
                               scalars->GetDataType(),
                               scalars->GetNumberOfTuples(),
                               scalars->GetNumberOfComponents(),
                               outputFormat);
}

//----------------------------------------------------------------------------
// Color type converters in anonymous namespace
namespace
{

#define vtkScalarsToColorsLuminance(r, g, b) \
    ((r)*0.30 + (g)*0.59 + (b)*0.11)

void vtkScalarsToColorsLuminanceToLuminance(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents)
{
  do
    {
    *outPtr++ = *inPtr;
    inPtr += numComponents;
    }
  while (--count);
}

void vtkScalarsToColorsLuminanceToRGB(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents)
{
  do
    {
    unsigned char l = *inPtr;
    outPtr[0] = l;
    outPtr[1] = l;
    outPtr[2] = l;
    inPtr += numComponents;
    outPtr += 3;
    }
  while (--count);
}

void vtkScalarsToColorsRGBToLuminance(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents)
{
  do
    {
    unsigned char r = inPtr[0];
    unsigned char g = inPtr[1];
    unsigned char b = inPtr[2];
    *outPtr++ = static_cast<unsigned char>(
                  vtkScalarsToColorsLuminance(r, g, b) + 0.5);
    inPtr += numComponents;
    }
  while (--count);
}

void vtkScalarsToColorsRGBToRGB(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents)
{
  do
    {
    outPtr[0] = inPtr[0];
    outPtr[1] = inPtr[1];
    outPtr[2] = inPtr[2];
    inPtr += numComponents;
    outPtr += 3;
    }
  while (--count);
}

void vtkScalarsToColorsLuminanceToLuminanceAlpha(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);

  do
    {
    outPtr[0] = inPtr[0];
    outPtr[1] = a;
    inPtr += numComponents;
    outPtr += 2;
    }
  while (--count);
}

void vtkScalarsToColorsLuminanceToRGBA(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);

  do
    {
    unsigned char l = inPtr[0];
    outPtr[0] = l;
    outPtr[1] = l;
    outPtr[2] = l;
    outPtr[3] = a;
    inPtr += numComponents;
    outPtr += 4;
    }
  while (--count);
}

void vtkScalarsToColorsRGBToLuminanceAlpha(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);

  do
    {
    unsigned char r = inPtr[0];
    unsigned char g = inPtr[1];
    unsigned char b = inPtr[2];
    outPtr[0] = static_cast<unsigned char>(
                  vtkScalarsToColorsLuminance(r, g, b) + 0.5);
    outPtr[1] = a;
    inPtr += numComponents;
    outPtr += 2;
    }
  while (--count);
}

void vtkScalarsToColorsRGBToRGBA(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);

  do
    {
    outPtr[0] = inPtr[0];
    outPtr[1] = inPtr[1];
    outPtr[2] = inPtr[2];
    outPtr[3] = a;
    inPtr += numComponents;
    outPtr += 4;
    }
  while (--count);
}


void vtkScalarsToColorsLuminanceAlphaToLuminanceAlpha(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  if (alpha >= 1)
    {
    do
      {
      outPtr[0] = inPtr[0];
      outPtr[1] = inPtr[1];
      inPtr += numComponents;
      outPtr += 2;
      }
    while (--count);
    }
  else
    {
    do
      {
      outPtr[0] = inPtr[0];
      outPtr[1] = static_cast<unsigned char>(inPtr[1]*alpha + 0.5);
      inPtr += numComponents;
      outPtr += 2;
      }
    while (--count);
    }
}

void vtkScalarsToColorsLuminanceAlphaToRGBA(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  if (alpha >= 1)
    {
    do
      {
      unsigned char l = inPtr[0];
      unsigned char a = inPtr[1];
      outPtr[0] = l;
      outPtr[1] = l;
      outPtr[2] = l;
      outPtr[3] = a;
      inPtr += numComponents;
      outPtr += 4;
      }
    while (--count);
    }
  else
    {
    do
      {
      unsigned char l = inPtr[0];
      unsigned char a = inPtr[1];
      outPtr[0] = l;
      outPtr[1] = l;
      outPtr[2] = l;
      outPtr[3] = static_cast<unsigned char>(a*alpha + 0.5);
      inPtr += numComponents;
      outPtr += 4;
      }
    while (--count);
    }
}

void vtkScalarsToColorsRGBAToLuminanceAlpha(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  do
    {
    unsigned char r = inPtr[0];
    unsigned char g = inPtr[1];
    unsigned char b = inPtr[2];
    unsigned char a = inPtr[3];
    outPtr[0] = static_cast<unsigned char>(
                  vtkScalarsToColorsLuminance(r, g, b) + 0.5);
    outPtr[1] = static_cast<unsigned char>(a*alpha + 0.5);
    inPtr += numComponents;
    outPtr += 2;
    }
  while (--count);
}

void vtkScalarsToColorsRGBAToRGBA(
  const unsigned char *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double alpha)
{
  if (alpha >= 1)
    {
    do
      {
      outPtr[0] = inPtr[0];
      outPtr[1] = inPtr[1];
      outPtr[2] = inPtr[2];
      outPtr[3] = inPtr[3];
      inPtr += numComponents;
      outPtr += 4;
      }
    while (--count);
    }
  else
    {
    do
      {
      outPtr[0] = inPtr[0];
      outPtr[1] = inPtr[1];
      outPtr[2] = inPtr[2];
      outPtr[3] = static_cast<unsigned char>(inPtr[3]*alpha + 0.5);
      inPtr += numComponents;
      outPtr += 4;
      }
    while (--count);
    }
}

//----------------------------------------------------------------------------

template<class T>
void vtkScalarsToColorsLuminanceToLuminance(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double l = inPtr[0];
    l = (l + shift)*scale;
    l = (l > minval ? l : minval);
    l = (l < maxval ? l : maxval);
    l += 0.5;
    outPtr[0] = static_cast<unsigned char>(l);
    inPtr += numComponents;
    outPtr += 1;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsLuminanceToRGB(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double l = inPtr[0];
    l = (l + shift)*scale;
    l = (l > minval ? l : minval);
    l = (l < maxval ? l : maxval);
    unsigned char lc = static_cast<unsigned char>(l + 0.5);
    outPtr[0] = lc;
    outPtr[1] = lc;
    outPtr[2] = lc;
    inPtr += numComponents;
    outPtr += 3;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsRGBToLuminance(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double r = inPtr[0];
    double g = inPtr[1];
    double b = inPtr[2];
    r = (r + shift)*scale;
    g = (g + shift)*scale;
    b = (b + shift)*scale;
    r = (r > minval ? r : minval);
    r = (r < maxval ? r : maxval);
    g = (g > minval ? g : minval);
    g = (g < maxval ? g : maxval);
    b = (b > minval ? b : minval);
    b = (b < maxval ? b : maxval);
    double l = vtkScalarsToColorsLuminance(r, g, b) + 0.5;
    outPtr[0] = static_cast<unsigned char>(l);
    inPtr += numComponents;
    outPtr += 1;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsRGBToRGB(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double r = inPtr[0];
    double g = inPtr[1];
    double b = inPtr[2];
    r = (r + shift)*scale;
    g = (g + shift)*scale;
    b = (b + shift)*scale;
    r = (r > minval ? r : minval);
    r = (r < maxval ? r : maxval);
    g = (g > minval ? g : minval);
    g = (g < maxval ? g : maxval);
    b = (b > minval ? b : minval);
    b = (b < maxval ? b : maxval);
    r += 0.5;
    g += 0.5;
    b += 0.5;
    outPtr[0] = static_cast<unsigned char>(r);
    outPtr[1] = static_cast<unsigned char>(g);
    outPtr[2] = static_cast<unsigned char>(b);
    inPtr += numComponents;
    outPtr += 3;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsLuminanceToLuminanceAlpha(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double l = inPtr[0];
    l = (l + shift)*scale;
    l = (l > minval ? l : minval);
    l = (l < maxval ? l : maxval);
    l += 0.5;
    outPtr[0] = static_cast<unsigned char>(l);
    outPtr[1] = a;
    inPtr += numComponents;
    outPtr += 2;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsLuminanceToRGBA(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double l = inPtr[0];
    l = (l + shift)*scale;
    l = (l > minval ? l : minval);
    l = (l < maxval ? l : maxval);
    unsigned char lc = static_cast<unsigned char>(l + 0.5);
    outPtr[0] = lc;
    outPtr[1] = lc;
    outPtr[2] = lc;
    outPtr[3] = a;
    inPtr += numComponents;
    outPtr += 4;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsRGBToLuminanceAlpha(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double r = inPtr[0];
    double g = inPtr[1];
    double b = inPtr[2];
    r = (r + shift)*scale;
    g = (g + shift)*scale;
    b = (b + shift)*scale;
    r = (r > minval ? r : minval);
    r = (r < maxval ? r : maxval);
    g = (g > minval ? g : minval);
    g = (g < maxval ? g : maxval);
    b = (b > minval ? b : minval);
    b = (b < maxval ? b : maxval);
    double l = vtkScalarsToColorsLuminance(r, g, b) + 0.5;
    outPtr[0] = static_cast<unsigned char>(l);
    outPtr[1] = a;
    inPtr += numComponents;
    outPtr += 2;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsRGBToRGBA(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double r = inPtr[0];
    double g = inPtr[1];
    double b = inPtr[2];
    r = (r + shift)*scale;
    g = (g + shift)*scale;
    b = (b + shift)*scale;
    r = (r > minval ? r : minval);
    r = (r < maxval ? r : maxval);
    g = (g > minval ? g : minval);
    g = (g < maxval ? g : maxval);
    b = (b > minval ? b : minval);
    b = (b < maxval ? b : maxval);
    r += 0.5;
    g += 0.5;
    b += 0.5;
    outPtr[0] = static_cast<unsigned char>(r);
    outPtr[1] = static_cast<unsigned char>(g);
    outPtr[2] = static_cast<unsigned char>(b);
    outPtr[3] = a;
    inPtr += numComponents;
    outPtr += 4;
    }
  while (--count);
}


template<class T>
void vtkScalarsToColorsLuminanceAlphaToLuminanceAlpha(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double l = inPtr[0];
    double a = inPtr[1];
    l = (l + shift)*scale;
    a = (a + shift)*scale;
    l = (l > minval ? l : minval);
    l = (l < maxval ? l : maxval);
    a = (a > minval ? a : minval);
    a = (a < maxval ? a : maxval);
    l += 0.5;
    a = a*alpha + 0.5;
    outPtr[0] = static_cast<unsigned char>(l);
    outPtr[1] = static_cast<unsigned char>(a);
    inPtr += numComponents;
    outPtr += 2;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsLuminanceAlphaToRGBA(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double l = inPtr[0];
    double a = inPtr[1];
    l = (l + shift)*scale;
    a = (a + shift)*scale;
    l = (l > minval ? l : minval);
    l = (l < maxval ? l : maxval);
    a = (a > minval ? a : minval);
    a = (a < maxval ? a : maxval);
    unsigned char lc = static_cast<unsigned char>(l + 0.5);
    a = a*alpha + 0.5;
    outPtr[0] = lc;
    outPtr[1] = lc;
    outPtr[2] = lc;
    outPtr[3] = static_cast<unsigned char>(a);
    inPtr += numComponents;
    outPtr += 4;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsRGBAToLuminanceAlpha(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double r = inPtr[0];
    double g = inPtr[1];
    double b = inPtr[2];
    double a = inPtr[3];
    r = (r + shift)*scale;
    g = (g + shift)*scale;
    b = (b + shift)*scale;
    a = (a + shift)*scale;
    r = (r > minval ? r : minval);
    r = (r < maxval ? r : maxval);
    g = (g > minval ? g : minval);
    g = (g < maxval ? g : maxval);
    b = (b > minval ? b : minval);
    b = (b < maxval ? b : maxval);
    a = (a > minval ? a : minval);
    a = (a < maxval ? a : maxval);
    a = a*alpha + 0.5;
    double l = vtkScalarsToColorsLuminance(r, g, b) + 0.5;
    outPtr[0] = static_cast<unsigned char>(l);
    outPtr[1] = static_cast<unsigned char>(a);
    inPtr += numComponents;
    outPtr += 2;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsRGBAToRGBA(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    double r = inPtr[0];
    double g = inPtr[1];
    double b = inPtr[2];
    double a = inPtr[3];
    r = (r + shift)*scale;
    g = (g + shift)*scale;
    b = (b + shift)*scale;
    a = (a + shift)*scale;
    r = (r > minval ? r : minval);
    r = (r < maxval ? r : maxval);
    g = (g > minval ? g : minval);
    g = (g < maxval ? g : maxval);
    b = (b > minval ? b : minval);
    b = (b < maxval ? b : maxval);
    a = (a > minval ? a : minval);
    a = (a < maxval ? a : maxval);
    r += 0.5;
    g += 0.5;
    b += 0.5;
    a = a*alpha + 0.5;
    outPtr[0] = static_cast<unsigned char>(r);
    outPtr[1] = static_cast<unsigned char>(g);
    outPtr[2] = static_cast<unsigned char>(b);
    outPtr[3] = static_cast<unsigned char>(a);
    inPtr += numComponents;
    outPtr += 4;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsMagnitudeToLuminance(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    int n = numComponents;
    double val = 0;
    do
      {
      double v = *inPtr++;
      val += v*v;
      }
    while (--n > 0);
    val = (sqrt(val) + shift)*scale;
    val = (val > minval ? val : minval);
    val = (val < maxval ? val : maxval);
    *outPtr++ = static_cast<unsigned char>(val + 0.5);
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsMagnitudeToLuminanceAlpha(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    int n = numComponents;
    double val = 0;
    do
      {
      double v = *inPtr++;
      val += v*v;
      }
    while (--n > 0);
    val = (sqrt(val) + shift)*scale;
    val = (val > minval ? val : minval);
    val = (val < maxval ? val : maxval);
    outPtr[0] = static_cast<unsigned char>(val + 0.5);
    outPtr[1] = a;
    outPtr += 2;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsMagnitudeToRGB(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale)
{
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    int n = numComponents;
    double val = 0;
    do
      {
      double v = *inPtr++;
      val += v*v;
      }
    while (--n > 0);
    val = (sqrt(val) + shift)*scale;
    val = (val > minval ? val : minval);
    val = (val < maxval ? val : maxval);
    unsigned char c = static_cast<unsigned char>(val + 0.5);
    outPtr[0] = c;
    outPtr[1] = c;
    outPtr[2] = c;
    outPtr += 3;
    }
  while (--count);
}

template<class T>
void vtkScalarsToColorsMagnitudeToRGBA(
  const T *inPtr, unsigned char *outPtr, vtkIdType count,
  int numComponents, double shift, double scale, double alpha)
{
  unsigned char a = static_cast<unsigned char>(alpha*255 + 0.5);
  static double minval = 0;
  static double maxval = 255.0;

  do
    {
    int n = numComponents;
    double val = 0;
    do
      {
      double v = *inPtr++;
      val += v*v;
      }
    while (--n > 0);
    val = (sqrt(val) + shift)*scale;
    val = (val > minval ? val : minval);
    val = (val < maxval ? val : maxval);
    unsigned char c = static_cast<unsigned char>(val + 0.5);
    outPtr[0] = c;
    outPtr[1] = c;
    outPtr[2] = c;
    outPtr[3] = a;
    outPtr += 4;
    }
  while (--count);
}

//----------------------------------------------------------------------------
unsigned char *vtkScalarsToColorsUnpackBits(void *inPtr, vtkIdType numValues)
{
  vtkIdType n = (numValues + 7) % 8;
  unsigned char *newPtr = new unsigned char [n];

  unsigned char *tmpPtr = newPtr;
  unsigned char *bitdata = static_cast<unsigned char *>(inPtr);
  for (vtkIdType i = 0; i < n; i += 8)
    {
    unsigned char b = *bitdata++;
    int j = 8;
    do
      {
      *tmpPtr++ = ((b >> (--j)) & 0x01);
      }
    while (j);
    }

  return newPtr;
}

// end anonymous namespace
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::MapColorsToColors(
  void *inPtr, unsigned char *outPtr, int inputDataType,
  int numberOfTuples, int numberOfComponents, int outputFormat)
{
  if (outputFormat < VTK_LUMINANCE || outputFormat > VTK_RGBA)
    {
    vtkErrorMacro(<< "MapScalarsToColors: unrecognized color format");
    return;
    }

  if (numberOfTuples <= 0)
    {
    return;
    }

  unsigned char *newPtr = 0;
  if (inputDataType == VTK_BIT)
    {
    newPtr = vtkScalarsToColorsUnpackBits(
      inPtr, numberOfTuples*numberOfComponents);
    inPtr = newPtr;
    inputDataType = VTK_UNSIGNED_CHAR;
    }

  double shift, scale;
  vtkScalarsToColorsComputeShiftScale(this, shift, scale);
  scale *= 255.0;

  double alpha = this->Alpha;
  if (alpha < 0) { alpha = 0; }
  if (alpha > 1) { alpha = 1; }

  if (inputDataType == VTK_UNSIGNED_CHAR &&
      static_cast<int>(shift*scale + 0.5) == 0 &&
      static_cast<int>((255 + shift)*scale + 0.5) == 255)
    {
    if (outputFormat == VTK_RGBA)
      {
      if (numberOfComponents == 1)
        {
        vtkScalarsToColorsLuminanceToRGBA(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      else if (numberOfComponents == 2)
        {
        vtkScalarsToColorsLuminanceAlphaToRGBA(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      else if (numberOfComponents == 3)
        {
        vtkScalarsToColorsRGBToRGBA(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      else
        {
        vtkScalarsToColorsRGBAToRGBA(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      }
    else if (outputFormat == VTK_RGB)
      {
      if (numberOfComponents < 3)
        {
        vtkScalarsToColorsLuminanceToRGB(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents);
        }
      else
        {
        vtkScalarsToColorsRGBToRGB(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents);
        }
      }
    else if (outputFormat == VTK_LUMINANCE_ALPHA)
      {
      if (numberOfComponents == 1)
        {
        vtkScalarsToColorsLuminanceToLuminanceAlpha(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      else if (numberOfComponents == 2)
        {
        vtkScalarsToColorsLuminanceAlphaToLuminanceAlpha(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      else if (numberOfComponents == 3)
        {
        vtkScalarsToColorsRGBToLuminanceAlpha(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      else
        {
        vtkScalarsToColorsRGBAToLuminanceAlpha(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents, alpha);
        }
      }
    else if (outputFormat == VTK_LUMINANCE)
      {
      if (numberOfComponents < 3)
        {
        vtkScalarsToColorsLuminanceToLuminance(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents);
        }
      else
        {
        vtkScalarsToColorsRGBToLuminance(
          static_cast<unsigned char*>(inPtr), outPtr,
          numberOfTuples, numberOfComponents);
        }
      }
    }
  else
    {
    // must apply shift scale and/or do type conversion
    if (outputFormat == VTK_RGBA)
      {
      if (numberOfComponents == 1)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsLuminanceToRGBA(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      else if (numberOfComponents == 2)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsLuminanceAlphaToRGBA(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      else if (numberOfComponents == 3)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsRGBToRGBA(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      else
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsRGBAToRGBA(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      }
    else if (outputFormat == VTK_RGB)
      {
      if (numberOfComponents < 3)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsLuminanceToRGB(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale));
          }
        }
      else
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsRGBToRGB(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale));
          }
        }
      }
    else if (outputFormat == VTK_LUMINANCE_ALPHA)
      {
      if (numberOfComponents == 1)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsLuminanceToLuminanceAlpha(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      else if (numberOfComponents == 2)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsLuminanceAlphaToLuminanceAlpha(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      else if (numberOfComponents == 3)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsRGBToLuminanceAlpha(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      else
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsRGBAToLuminanceAlpha(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale, alpha));
          }
        }
      }
    else if (outputFormat == VTK_LUMINANCE)
      {
      if (numberOfComponents < 3)
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsLuminanceToLuminance(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale));
          }
        }
      else
        {
        switch (inputDataType)
          {
          vtkTemplateAliasMacro(
            vtkScalarsToColorsRGBToLuminance(
              static_cast<VTK_TT*>(inPtr), outPtr,
              numberOfTuples, numberOfComponents, shift, scale));
          }
        }
      }
    }

  if (newPtr)
    {
    delete [] newPtr;
    }
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::MapScalarsThroughTable2(
  void *inPtr, unsigned char *outPtr, int inputDataType,
  int numberOfTuples, int numberOfComponents, int outputFormat)
{
  if (outputFormat < VTK_LUMINANCE || outputFormat > VTK_RGBA)
    {
    vtkErrorMacro(<< "MapScalarsThroughTable2: unrecognized color format");
    return;
    }

  if (numberOfTuples <= 0)
    {
    return;
    }

  unsigned char *newPtr = 0;
  if (inputDataType == VTK_BIT)
    {
    newPtr = vtkScalarsToColorsUnpackBits(
      inPtr, numberOfTuples*numberOfComponents);
    inPtr = newPtr;
    inputDataType = VTK_UNSIGNED_CHAR;
    }

  double shift, scale;
  vtkScalarsToColorsComputeShiftScale(this, shift, scale);
  scale *= 255.0;

  double alpha = this->Alpha;
  if (alpha < 0) { alpha = 0; }
  if (alpha > 1) { alpha = 1; }

  if (this->UseMagnitude && numberOfComponents > 1)
    {
    if (outputFormat == VTK_RGBA)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsMagnitudeToRGBA(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale, alpha));
        }
      }
    else if (outputFormat == VTK_RGB)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsMagnitudeToRGB(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale));
        }
      }
    else if (outputFormat == VTK_LUMINANCE_ALPHA)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsMagnitudeToLuminanceAlpha(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale, alpha));
        }
      }
    else if (outputFormat == VTK_LUMINANCE)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsMagnitudeToLuminance(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale));
        }
      }
    }
  else if (inputDataType == VTK_UNSIGNED_CHAR &&
      static_cast<int>(shift*scale + 0.5) == 0 &&
      static_cast<int>((255 + shift)*scale + 0.5) == 255)
    {
    if (outputFormat == VTK_RGBA)
      {
      vtkScalarsToColorsLuminanceToRGBA(
        static_cast<unsigned char*>(inPtr), outPtr,
        numberOfTuples, numberOfComponents, alpha);
      }
    else if (outputFormat == VTK_RGB)
      {
      vtkScalarsToColorsLuminanceToRGB(
        static_cast<unsigned char*>(inPtr), outPtr,
        numberOfTuples, numberOfComponents);
      }
    else if (outputFormat == VTK_LUMINANCE_ALPHA)
      {
      vtkScalarsToColorsLuminanceToLuminanceAlpha(
        static_cast<unsigned char*>(inPtr), outPtr,
        numberOfTuples, numberOfComponents, alpha);
      }
    else if (outputFormat == VTK_LUMINANCE)
      {
      vtkScalarsToColorsLuminanceToLuminance(
        static_cast<unsigned char*>(inPtr), outPtr,
        numberOfTuples, numberOfComponents);
      }
    }
  else
    {
    // must apply shift scale and/or do type conversion
    if (outputFormat == VTK_RGBA)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsLuminanceToRGBA(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale, alpha));
        }
      }
    else if (outputFormat == VTK_RGB)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsLuminanceToRGB(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale));
        }
      }
    else if (outputFormat == VTK_LUMINANCE_ALPHA)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsLuminanceToLuminanceAlpha(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale, alpha));
        }
      }
    else if (outputFormat == VTK_LUMINANCE)
      {
      switch (inputDataType)
        {
        vtkTemplateAliasMacro(
          vtkScalarsToColorsLuminanceToLuminance(
            static_cast<VTK_TT*>(inPtr), outPtr,
            numberOfTuples, numberOfComponents, shift, scale));
        }
      }
    }

  if (newPtr)
    {
    delete [] newPtr;
    }
}

//----------------------------------------------------------------------------
vtkUnsignedCharArray *vtkScalarsToColors::ConvertUnsignedCharToRGBA(
  vtkUnsignedCharArray *colors, int numComp, int numTuples)
{
  if ( numComp == 4 && this->Alpha >= 1.0 )
    {
    colors->Register(this);
    return colors;
    }

  unsigned char *cptr = colors->GetPointer(0);
  vtkUnsignedCharArray *newColors = vtkUnsignedCharArray::New();
  newColors->SetNumberOfComponents(4);
  newColors->SetNumberOfTuples(numTuples);
  unsigned char *nptr = newColors->GetPointer(0);
  double alpha = this->Alpha;
  alpha = (alpha > 0 ? alpha : 0);
  alpha = (alpha < 1 ? alpha : 1);

  if (numTuples <= 0)
    {
    return newColors;
    }

  switch (numComp)
    {
    case 1:
      vtkScalarsToColorsLuminanceToRGBA(
        cptr, nptr, numTuples, numComp, alpha);
      break;

    case 2:
      vtkScalarsToColorsLuminanceAlphaToRGBA(
        cptr, nptr, numTuples, numComp, alpha);
      break;

    case 3:
      vtkScalarsToColorsRGBToRGBA(
        cptr, nptr, numTuples, numComp, alpha);
      break;

    case 4:
      vtkScalarsToColorsRGBAToRGBA(
        cptr, nptr, numTuples, numComp, alpha);
      break;

    default:
      vtkErrorMacro(<<"Cannot convert colors");
      return NULL;
    }
  
  return newColors;
}

//----------------------------------------------------------------------------
void vtkScalarsToColors::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Alpha: " << this->Alpha << endl;
  if (this->VectorMode == vtkScalarsToColors::MAGNITUDE)
    {
    os << indent << "VectorMode: Magnitude\n";
    }
  else if (this->VectorMode == vtkScalarsToColors::RGBCOLORS)
    {
    os << indent << "VectorMode: RGBColors\n";
    }
  else
    {
    os << indent << "VectorMode: Component\n";
    os << indent << "VectorComponent: " << this->VectorComponent << endl;
    }
}
