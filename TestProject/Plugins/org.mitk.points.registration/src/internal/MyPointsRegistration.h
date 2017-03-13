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


#ifndef MyPointsRegistration_h
#define MyPointsRegistration_h

#include <berryISelectionListener.h>

#include <QmitkAbstractView.h>

#include "ui_MyPointsRegistrationControls.h"
#include <mitkPointSet.h>
#include <mitkStandaloneDataStorage.h>
#include <itkImage.h>
#include <QWidget>
#include <mitkImage.h>
#include <eigen3/Eigen/Core>
#include <QTableWidget>

/**
  \brief MyPointsRegistration

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class MyPointsRegistration : public QmitkAbstractView
{
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

  public:

    static const std::string VIEW_ID;
   
  protected slots:

    /// \brief Called when the user clicks the GUI button
    void CreatePointSet1();
    void createPointSet2();
    
    void performRegistration();
    
	Eigen::Matrix<double, 4, 4> point_based_registration ();
    Eigen::Vector3d ComputeSetPointCentroid(mitk::PointSet::Pointer, int);
    Eigen::MatrixXd PointSetDeviated(mitk::PointSet::Pointer, int ,Eigen::Vector3d );
    double checkFiducialError();
    

  protected:
    
    mitk::PointSet::Pointer m_fixedPointSet;
    mitk::PointSet::Pointer m_movingPointSet; 
    Eigen::MatrixXd registeredSet_matrix;
    Eigen::MatrixXd  transformation_matrix;
    virtual void CreateQtPartControl(QWidget *parent) override;

    virtual void SetFocus() override;

    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged( berry::IWorkbenchPart::Pointer source,
                                     const QList<mitk::DataNode::Pointer>& nodes ) override;

    Ui::MyPointsRegistrationControls m_Controls;

};

#endif // MyPointsRegistration_h
