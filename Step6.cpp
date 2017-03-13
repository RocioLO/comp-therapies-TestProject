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
#include "Step6.h"
#include "QmitkRenderWindow.h"
#include "QmitkSliceWidget.h"
#include "mitkProperties.h"
#include "mitkRenderingManager.h"
#include "mitkPointSet.h"
#include "mitkPointSetDataInteractor.h"
#include "mitkImageAccessByItk.h"
#include "mitkRenderingManager.h"
#include <mitkIOUtil.h>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QPushButton>
#include <QVBoxLayout>

#include <vtkImageData.h>
#include <vtkMarchingCubes.h>
#include <vtkSTLWriter.h>

#include <QmitkStdMultiWidget.h>
		 

//##Documentation
//## @brief Start region-grower at interactively added points
Step6::Step6(int argc, char *argv[], QWidget *parent) 
: QWidget(parent)
{
  // load data as in the previous steps; a reference to the first loaded
  // image is kept in the member m_FirstImage and used as input for the
  // region growing
  Load(argc, argv);
}
void Step6::Initialize()
{
  // setup the widgets as in the previous steps, but with an additional
  // QVBox for a button to start the segmentation
  this->SetupWidgets();
  // Create controlsParent widget with horizontal layout
  QWidget *controlsParent = new QWidget(this);
  this->layout()->addWidget(controlsParent);
  QHBoxLayout *hlayout = new QHBoxLayout(controlsParent);
  hlayout->setSpacing(2);
  QLabel *labelThresholdMin = new QLabel("Lower Threshold:", controlsParent);
  hlayout->addWidget(labelThresholdMin);
  m_LineEditThresholdMin = new QLineEdit("-1000", controlsParent);
  hlayout->addWidget(m_LineEditThresholdMin);
  QLabel *labelThresholdMax = new QLabel("Upper Threshold:", controlsParent);
  hlayout->addWidget(labelThresholdMax);
  m_LineEditThresholdMax = new QLineEdit("-400", controlsParent);
  hlayout->addWidget(m_LineEditThresholdMax);
  // create button to start the segmentation and connect its clicked()
  // signal to method StartRegionGrowing
  QPushButton *startButton = new QPushButton("start region growing", controlsParent);
  hlayout->addWidget(startButton);
  
  connect(startButton, SIGNAL(clicked()), this, SLOT(StartRegionGrowing()));
  if (m_FirstImage.IsNull())
    startButton->setEnabled(false);
  // as in Step5, create PointSet (now as a member m_Seeds) and
  // associate a interactor to it
  m_Seeds = mitk::PointSet::New();
  mitk::DataNode::Pointer pointSetNode = mitk::DataNode::New();
  pointSetNode->SetData(m_Seeds);
  pointSetNode->SetProperty("layer", mitk::IntProperty::New(2));
  m_DataStorage->Add(pointSetNode);
  // Create PointSetDataInteractor
  mitk::PointSetDataInteractor::Pointer interactor = mitk::PointSetDataInteractor::New();
  interactor->LoadStateMachine("PointSet.xml");
  interactor->SetEventConfig("PointSetConfig.xml");
  interactor->SetDataNode(pointSetNode);
}
int Step6::GetThresholdMin()
{
  return m_LineEditThresholdMin->text().toInt();
}
int Step6::GetThresholdMax()
{
  return m_LineEditThresholdMax->text().toInt();
}
void Step6::StartRegionGrowing()
{
  AccessByItk_1(m_FirstImage, RegionGrowing, this);
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  
  std::cout << "7";
  if (m_ResultImage.IsNotNull())
  {
    m_ResultNode->SetProperty("volumerendering", mitk::BoolProperty::New(false));
    vtkMarchingCubes *surfaceCreator = vtkMarchingCubes::New();
    surfaceCreator->SetInputData(m_ResultImage->GetVtkImageData());
    surfaceCreator->SetValue(0, 1);
    surfaceCreator->Update();
    mitk::Surface::Pointer surface = mitk::Surface::New();
    surface->SetVtkPolyData(surfaceCreator->GetOutput()); // VTK6_TODO
    mitk::DataNode::Pointer surfaceNode = mitk::DataNode::New();
    surfaceNode->SetData(surface);
    m_DataStorage->Add(surfaceNode);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    std::cout << "8";
    surfaceCreator->Delete();
  }
  std::cout << "9";
}
void Step6::Load(int argc, char *argv[])
{
  //*************************************************************************
  // Part I: Basic initialization
  //*************************************************************************
  m_DataStorage = mitk::StandaloneDataStorage::New();
  //*************************************************************************
  // Part II: Create some data by reading files
  //*************************************************************************
  int i;
  for (i = 1; i < argc; ++i)
  {
    // For testing
    if (strcmp(argv[i], "-testing") == 0)
      continue;
    // Load datanode (eg. many image formats, surface formats, etc.)
    mitk::StandaloneDataStorage::SetOfObjects::Pointer dataNodes = mitk::IOUtil::Load(argv[i], *m_DataStorage);
    if (dataNodes->empty())
    {
      fprintf(stderr, "Could not open file %s \n\n", argv[i]);
      exit(2);
    }
    mitk::Image::Pointer image = dynamic_cast<mitk::Image *>(dataNodes->at(0)->GetData());
    if ((m_FirstImage.IsNull()) && (image.IsNotNull()))
      m_FirstImage = image;
  }
}
void Step6::SetupWidgets()
{
  //*************************************************************************
  // Part I: Create windows and pass the datastorage to it
  //*************************************************************************
  // Create toplevel widget with vertical layout
  QVBoxLayout *vlayout = new QVBoxLayout(this);
  vlayout->setMargin(0);
  //vlayout->setSpacing(2);
  // Create viewParent widget with horizontal layout
  QWidget *viewParent = new QWidget(this);
  vlayout->addWidget(viewParent);
  QHBoxLayout *hlayout = new QHBoxLayout(viewParent);
  hlayout->setMargin(0);
  //hlayout->setSpacing(2);
  //*************************************************************************
  // Step 8: Using QmitkStdMultiWidget creation and initialization
  //*************************************************************************
  QmitkStdMultiWidget *multiWidget = new QmitkStdMultiWidget(viewParent);
  hlayout->addWidget(multiWidget);
  // Tell the multiWidget which DataStorage to render
  multiWidget->SetDataStorage(m_DataStorage);
  // Initialize views as axial, sagittal, coronar (from
  // top-left to bottom)
  mitk::TimeGeometry::Pointer geo = m_DataStorage->ComputeBoundingGeometry3D(m_DataStorage->GetAll());
  mitk::RenderingManager::GetInstance()->InitializeViews(geo);
  // Initialize bottom-right view as 3D view
  multiWidget->GetRenderWindow4()->GetRenderer()->SetMapperID(mitk::BaseRenderer::Standard3D);
  // Enable standard handler for levelwindow-slider
  multiWidget->EnableStandardLevelWindow();
  // Add the displayed views to the DataStorage to see their positions in 2D and 3D
  multiWidget->AddDisplayPlaneSubTree();
  multiWidget->AddPlanesToDataStorage();
  multiWidget->SetWidgetPlanesVisibility(true);
  
}
