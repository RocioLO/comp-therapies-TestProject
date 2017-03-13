/*===================================================================
The Medical Imaging Interaction Toolkit (MITK)
Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.
See LICENSE.txt or http://www.mitk.org for details.
===================================================================*/
#include "Step6RegionGrowing.txx"
#include <mitkInstantiateAccessFunctions.h>
#define InstantiateAccessFunction_RegionGrowing(pixelType, dim)                                                        \
  \
template void                                                                                                          \
    RegionGrowing(itk::Image<pixelType, dim> *, Step6 *);
InstantiateAccessFunctionForFixedDimension(RegionGrowing, 3)
  
