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


// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk

#include "QmitkRegisterClasses.h"
#include "QmitkRenderWindow.h"
#include "QmitkSliceWidget.h"
#include <mitkDataNodeSelection.h>
// Qt
#include <QMessageBox>
#include <QApplication>
#include <QHBoxLayout>

//mitk image
#include <mitkImage.h>
#include "mitkNodePredicateDataType.h"
#include "mitkProperties.h"
#include "mitkRenderingManager.h"
#include "mitkStandaloneDataStorage.h"
#include "mitkPointSet.h"
#include "mitkPointOperation.h"
#include "mitkInteractionConst.h"
#include "mitkPointSetDataInteractor.h"
#include <mitkIOUtil.h>

#include <iostream>
#include <iomanip>
#include <string>

#include <itksys/SystemTools.hxx>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Dense>
#include <math.h>
#include <vtkMath.h>
#include <QmitkPointListWidget.h>

#include "MyPointsRegistration.h"
using namespace Eigen;

using namespace std;


const std::string MyPointsRegistration::VIEW_ID = "org.mitk.views.mypointsregistration";

void MyPointsRegistration::SetFocus()
{
  m_Controls.buttonPointSet1->setFocus();
}

void MyPointsRegistration::CreateQtPartControl( QWidget *parent )
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi( parent );
  
  connect( m_Controls.buttonPointSet1, SIGNAL(clicked()), this, SLOT(CreatePointSet1()) );
  connect( m_Controls.buttonPointSet2, SIGNAL(clicked()), this, SLOT(createPointSet2()) );
  connect( m_Controls.registrationButton, SIGNAL(clicked()), this, SLOT(performRegistration()) );
}

void MyPointsRegistration::OnSelectionChanged( berry::IWorkbenchPart::Pointer /*source*/,
                                             const QList<mitk::DataNode::Pointer>& nodes )
{
  // iterate all selected objects, adjust warning visibility
  foreach( mitk::DataNode::Pointer node, nodes )
  {
    if( node.IsNotNull() && dynamic_cast<mitk::Image*>(node->GetData()) )
    {
      m_Controls.labelWarning->setVisible( false );
      m_Controls.buttonPointSet1->setEnabled( true );
      return;
    }
  }

  m_Controls.labelWarning->setVisible( true );
  m_Controls.buttonPointSet1->setEnabled( true );
  m_Controls.buttonPointSet2->setEnabled( true );
 }


void MyPointsRegistration::CreatePointSet1()
{
  // Point Set 1 will be the fixed set of points to which we will "adjust" Point Set 2
  //We create an empty pointset
  m_fixedPointSet = mitk::PointSet::New();
  
  //We create a new data tree node
  mitk::DataNode::Pointer pointSetNode = mitk::DataNode::New();
  
  // Store the point set in the DataNode and edit its properties
  pointSetNode->SetData(m_fixedPointSet); 
  pointSetNode->SetProperty("name", mitk::StringProperty::New("PointSet1"));
  pointSetNode->SetColor( 0.0, 1.0, 0.0 );
  
  // Add the node to the tree
  //datastorage->Add(pointSetNode);
  this->GetDataStorage()->Add(pointSetNode);
  // Create PointSetDataInteractor
  
  mitk::PointSetDataInteractor::Pointer interactor = mitk::PointSetDataInteractor::New();
  // Set the StateMachine pattern that describes the flow of the interactions
  interactor->LoadStateMachine("PointSet.xml");
  
  // Set the configuration file, which describes the user interactions that trigger actions
  // in this file SHIFT + LeftClick triggers add Point, but by modifying this file,
  // it could as well be changes to any other user interaction.
  interactor->SetEventConfig("PointSetConfig.xml");
  
  interactor->SetDataNode(pointSetNode); // Assign the pointSetNode to the interactor 
  
}
 
 
void MyPointsRegistration::createPointSet2 ()
{
	//In order to create only one set of points with name Point Set1 we set it to false once used
  m_Controls.buttonPointSet1->setEnabled(false);
  // Point Set 2 will be the moving set of points.
  
  //We create an empty pointset
  m_movingPointSet = mitk::PointSet::New();
  
  //We create a new data tree node
  mitk::DataNode::Pointer pointSetNode2 = mitk::DataNode::New();
  
  // Store the point set in the DataNode and edit its properties
  pointSetNode2->SetData(m_movingPointSet); 
  pointSetNode2->SetProperty("name", mitk::StringProperty::New("PointSet2"));
  pointSetNode2->SetColor( 0.0, 0.0, 1.0 );
  
  // Add the node to the tree
  //datastorage->Add(pointSetNode);
  this->GetDataStorage()->Add(pointSetNode2);
  // Create PointSetDataInteractor
  
  mitk::PointSetDataInteractor::Pointer interactor2 = mitk::PointSetDataInteractor::New();
  // Set the StateMachine pattern that describes the flow of the interactions
  interactor2->LoadStateMachine("PointSet.xml");
  // Set the configuration file, which describes the user interactions that trigger actions
  // in this file SHIFT + LeftClick triggers add Point, but by modifying this file,
  // it could as well be changes to any other user interaction.
  interactor2->SetEventConfig("PointSetConfig.xml");
  
  interactor2->SetDataNode(pointSetNode2); // Assign the pointSetNode to the interactor,
										   
  
}


void MyPointsRegistration::performRegistration()
{
	m_Controls.buttonPointSet2->setEnabled(false);
	//Registration from Point Set 2 (moving) to Point Set 1 (fixed)
	if (m_movingPointSet.IsNotNull() && m_fixedPointSet.IsNotNull()&& m_movingPointSet->GetSize() != 0 && m_fixedPointSet->GetSize() != 0)
    {
		Matrix<double, 4, 4> matrix_transformation; matrix_transformation << 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0;
        
        /*Matrix<double, 4, 4> registeredSet_matrix; registeredSet_matrix << 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0,
																			 0.0, 0.0, 0.0, 0.0;*/
        
		matrix_transformation = point_based_registration();
		//Print transformation in console to check
		//printf("Transformation Matrix computed sucessfully! \n");
		/*printf("Transformation Matrix is: %f, %f, %f , %f\n, %f, %f, %f, %f \n %f, %f, %f, %f \n, %f, %f, %f , %f\n", 
		matrix_transformation(0,0), matrix_transformation(0,1),matrix_transformation(0,2), matrix_transformation(0,3), 
		matrix_transformation(1,0), matrix_transformation(1,1), matrix_transformation(1,2), matrix_transformation(1,3), 
		matrix_transformation(2,0), matrix_transformation(2,1), matrix_transformation(2,2), matrix_transformation(2,3),
		matrix_transformation(3,0), matrix_transformation(3,1), matrix_transformation(3,2), matrix_transformation(3,3));*/

		//Point Set 2 registered to Point Set 1 
		mitk::Point3D point;  point[0]=0; point[1]=0; point[2]=0; 
		
		Matrix<double, 3, Dynamic> PointSet2_matrix (3, m_movingPointSet->GetSize());
		
		//This will convert PointSet2 in a matrix so that do calculations faster
		for (int pointInSet=0; pointInSet<m_movingPointSet->GetSize(); ++pointInSet)
		{
			
			point = m_movingPointSet->GetPoint(pointInSet);
			
			PointSet2_matrix(0,pointInSet) = point[0];
			PointSet2_matrix(1,pointInSet) = point[1];
			PointSet2_matrix(2,pointInSet) = point[2];
		}
		/*printf(" PointSet2 matrix is: %f, %f, %f, %f \n, %f, %f, %f, %f \n %f, %f, %f, %f \n ", 
		PointSet2_matrix(0,0), PointSet2_matrix(0,1),PointSet2_matrix(0,2), PointSet2_matrix(0,3),
		PointSet2_matrix(1,0), PointSet2_matrix(1,1), PointSet2_matrix(1,2), PointSet2_matrix(1,3),
		PointSet2_matrix(2,0), PointSet2_matrix(2,1), PointSet2_matrix(2,2), PointSet2_matrix(2,3));*/

		//We need squared matrices, so I create a temporal auxiliar expanding PointSet2_matrix row:
		Matrix<double, 4, Dynamic> aux_matrix(4, m_movingPointSet->GetSize());
		
		aux_matrix << PointSet2_matrix,
			          MatrixXd::Ones(1, m_movingPointSet->GetSize());
		
		
		registeredSet_matrix = matrix_transformation*aux_matrix;
		
		//We return the points from eigen matrix to Point Set
		mitk::Point3D point_transf;  point_transf[0]=0; point_transf[1]=0; point_transf[2]=0; 
		for (int pointInSet=0; pointInSet<  m_movingPointSet->GetSize(); ++pointInSet)
		{
			point_transf [0] = registeredSet_matrix(0,pointInSet);
			point_transf [1] = registeredSet_matrix(1,pointInSet); 
			point_transf [2] = registeredSet_matrix(2,pointInSet);
			
			m_movingPointSet->InsertPoint(pointInSet,point_transf); //Se puede sobreescribir un PointSet?
																 // Actualizar data-storage?
			//printf("Nuevos puntos del set son: %f, %f, %f \n", point_transf[0], point_transf[1], point_transf[2]);													 
		}	
		m_movingPointSet->GetTimeGeometry()->Update();
		m_movingPointSet->UpdateOutputInformation();
		m_movingPointSet->OnPointSetChange();
		m_movingPointSet->Modified();
		
		m_Controls.tableWidget->setRowCount(4);
		m_Controls.tableWidget->setColumnCount(4);
		m_Controls.tableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
		m_Controls.tableWidget->setShowGrid(true);
		
		//Show results in Qt Table
		for (int i=0; i<4; ++i)
		{
			for (int j=0; j<4; ++j)
			{
				double value = matrix_transformation(i,j);
				
				m_Controls.tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(value)));  
			}
		}
		
		//We check the registration error:
		double fiducialError = checkFiducialError();
		//printf("El error de registro RMSE es: %f", fiducialError);
		m_Controls.label_3->setText(QString::number(fiducialError));
				
	}
	else
	{
		printf("Check that the Point Sets are not empty");
	}
	
}

//Eigen::MatrixXd MyPointsRegistration::point_based_registration()
Eigen::Matrix<double, 4, 4> MyPointsRegistration::point_based_registration()
{   
	//Local variables initialization
	Vector3d moving_centroid(0.000, 0.000, 0.000); 
	Vector3d fixed_centroid(0.000, 0.000, 0.000); 
	
	int NumberOfFixedPoints = 0;
	int NumberOfMovingPoints = 0;
	
	/*Matrix<double, 4, 4> matrix_transformation; matrix_transformation << 0.0, 0.0, 0.0, 0.0,
																		 0.0, 0.0, 0.0, 0.0,
																		 0.0, 0.0, 0.0, 0.0,
																		 0.0, 0.0, 0.0, 0.0;
        */
	
	// We check if the number of points is equal in both sets, if not registration can't be performed.	
	NumberOfFixedPoints = m_fixedPointSet->GetSize();
	NumberOfMovingPoints = m_movingPointSet->GetSize();
	//printf ("Size Fixed Set: %d \n", NumberOfFixedPoints);
	//printf ("Size Moving Set: %d \n", NumberOfMovingPoints);
	
	if(NumberOfFixedPoints != NumberOfMovingPoints)
	{
		printf("Error: different number of points between fixed and moving set of points");
		
	}
	
	// Step 1: Deviations from centroid
	moving_centroid = ComputeSetPointCentroid(m_movingPointSet, NumberOfMovingPoints);
	
	fixed_centroid = ComputeSetPointCentroid(m_fixedPointSet, NumberOfFixedPoints);
	
	
	//printf("The fixed centroid is: %f, %f, %f \n",fixed_centroid(0), fixed_centroid(1), fixed_centroid(2) );
	//printf("The moving centroid is: %f, %f, %f \n",moving_centroid(0), moving_centroid(1), fixed_centroid(2) );
	
	MatrixXd q_fixed = PointSetDeviated(m_fixedPointSet, NumberOfFixedPoints, fixed_centroid);
	MatrixXd q_moving = PointSetDeviated(m_movingPointSet, NumberOfMovingPoints, moving_centroid);
	//printf("Matrix q_fixed is: %f, %f, %f \n, %f, %f, %f \n %f, %f, %f \n ", q_fixed(0,0), q_fixed(0,1),q_fixed(0,2), q_fixed(1,0), q_fixed(1,1), q_fixed(1,2), q_fixed(2,0), q_fixed(2,1), q_fixed(2,2));

	//Step 2: Calculate the 3x3 matrix
	Matrix<double, 3, 3> H_matrix = q_moving*q_fixed.transpose();
	//printf("Matrix H is: %f, %f, %f \n, %f, %f, %f \n %f, %f, %f \n ", H_matrix(0,0), H_matrix(0,1),H_matrix(0,2), H_matrix(1,0), H_matrix(1,1), H_matrix(1,2), H_matrix(2,0), H_matrix(2,1), H_matrix(2,2));
	//Step 3: Singular Value Decomposition of H
	JacobiSVD<MatrixXd> svd(H_matrix, ComputeThinU | ComputeThinV);
	Matrix<double, 3, 3> U_matrix = svd.matrixU();
	Matrix<double, 3, 3> V_matrix = svd.matrixV();
	
	//Step 4: R = V*U(transposed)
	Matrix<double, 3, 3> Rotation_matrix = V_matrix*U_matrix.transpose();
	//printf("Matrix Rotation is: %f, %f, %f \n, %f, %f, %f \n %f, %f, %f \n ", Rotation_matrix(0,0), Rotation_matrix(0,1),Rotation_matrix(0,2), Rotation_matrix(1,0), Rotation_matrix(1,1), Rotation_matrix(1,2), Rotation_matrix(2,0), Rotation_matrix(2,1), Rotation_matrix(2,2));

	// Check if the determinant of Rotation_matrix is +1 or -1:
	if (Rotation_matrix.determinant() == -1)
	{
		printf ("Error while computing Rotation Matrix");
	}
		
	
	//Step 5: 
	//Vector3d Transaltion_matrix = fixed_centroid.transpose() - Rotation_matrix*moving_centroid;
	Vector3d vector_aux =  Rotation_matrix  * moving_centroid;
	Vector3d Transaltion_matrix = fixed_centroid - vector_aux;
	//printf("El vector de translacion es: %f, %f, %f \n", Transaltion_matrix(0), Transaltion_matrix(1), Transaltion_matrix(2) );
	//Step 6: Transformation matrix will be a 4x4 matrix concatenating rotation and translation and adding
	// a row of: 0 0 0 1
	
	Matrix<double, 4, 4> transformation_matrix;
	transformation_matrix << Rotation_matrix, 
							 Transaltion_matrix,
							 MatrixXd::Zero(1,3),
							 MatrixXd::Identity(1,1);
							 
	return transformation_matrix;
	
}


Eigen::Vector3d MyPointsRegistration::ComputeSetPointCentroid(mitk::PointSet::Pointer pointSet, int numberOfPoints)
{
	mitk::Point3D point; point[0]=0; point[1]=0; point[2]=0; 
	double x =0.0;
	double y =0.0;
	double z =0.0;
	Vector3d centroid(0,0,0);
	
	for (int pointInSet=0; pointInSet< numberOfPoints; ++pointInSet)
	{
		point = pointSet->GetPoint(pointInSet);
		x += point[0];
		y += point[1];
		z += point[2];
		//printf("point son: %f, %f, %f \n", point[0], point[1], point[2]);
	}
	//printf("la suma de los x es: %f \n", x);
	
	centroid(0) = x/numberOfPoints; 
	centroid(1) = y/numberOfPoints;
	centroid(2) = z/numberOfPoints;
	
	return centroid;
}


Eigen::MatrixXd MyPointsRegistration::PointSetDeviated(mitk::PointSet::Pointer pointSet, int numberOfPoints, Vector3d centroid)
{
	//Every column in the matrix is a different 3D Point
	mitk::Point3D pointDeviated; pointDeviated[0]=0; pointDeviated[1]=0; pointDeviated[2]=0; 
	mitk::Point3D point;  point[0]=0; point[1]=0; point[2]=0; 
	Matrix<double, 3, Dynamic> fixed_deviated_matrix_points (3, numberOfPoints);

	for (int pointInSet=0; pointInSet< numberOfPoints; ++pointInSet)
	{
		//We center the set of points in the origin
	    point = pointSet->GetPoint(pointInSet);
		pointDeviated[0] = point[0] - centroid[0];
		pointDeviated[1] = point[1] - centroid[1];
		pointDeviated[2] = point[2] - centroid[2];
		//printf("Puntos cenrados: %f, %f, %f \n",pointDeviated[0], pointDeviated[1], pointDeviated[2] );
		
		fixed_deviated_matrix_points(0,pointInSet) = pointDeviated[0];//insert estaba en la documentacion de sparseMatrix, para
		fixed_deviated_matrix_points(1,pointInSet) = pointDeviated[1];//matrices cualquiera no estoy segura si es así o con =
		fixed_deviated_matrix_points(2,pointInSet) = pointDeviated[2];
	}
	//printf("Set de puntos desviados en matriz: %f, %f, %f \n", fixed_deviated_matrix_points(0,0), fixed_deviated_matrix_points(1,0), fixed_deviated_matrix_points(2,0));

	return fixed_deviated_matrix_points;
}


double MyPointsRegistration::checkFiducialError()
{
	mitk::Point3D point;  point[0]=0; point[1]=0; point[2]=0; 
	Matrix<double, 3, Dynamic> PointSet1_matrix (3, m_movingPointSet->GetSize());
	
	//This will convert PointSet1 in a matrix so that do calculations faster
	for (int pointInSet=0; pointInSet<m_fixedPointSet->GetSize(); ++pointInSet)
	{
		//We center the set of points in the origin
		point = m_fixedPointSet->GetPoint(pointInSet);
		
		PointSet1_matrix(0,pointInSet) = point[0];
		PointSet1_matrix(1,pointInSet) = point[1];
		PointSet1_matrix(2,pointInSet) = point[2];
	}
	//We need both matrices to be same size, as the registered set is a 4x4 matrix we 
	//convert point set 1 matrix into a 4x4 as well
	
	Matrix<double, 4, Dynamic> PointSet1_4x4matrix(4, m_fixedPointSet->GetSize());
	
	PointSet1_4x4matrix << PointSet1_matrix,
						   MatrixXd::Ones(1, m_fixedPointSet->GetSize());
	
    Array<double, 1, Dynamic > error_mm(1, m_fixedPointSet->GetSize());
	Vector4d aux_norm(0.000, 0.000, 0.000, 0.000);
	
	
	for (int i=0; i<m_fixedPointSet->GetSize(); ++i)
	{
		
		for(int j=0; j<4 ; ++j)
		{
			
			aux_norm(j) = PointSet1_4x4matrix(j,i) - registeredSet_matrix(j,i);
			
		    
		}
		error_mm(i) =  aux_norm.norm();
	}
	printf("vector con normales es %f, %f, %f, %f \n",error_mm(0),error_mm(1), error_mm(2),error_mm(3) );
	
	double count = 0.000;
	for (int i=0; i<m_fixedPointSet->GetSize(); ++i)
	{
		count += error_mm(i);
	}
	
	double error_rms_mm = count/m_fixedPointSet->GetSize();
	
	return sqrt(error_rms_mm);
}

